import csv
import os
import shutil
import subprocess
import sys
from pathlib import PureWindowsPath
import numpy as np
from rdkit import Chem
from rdkit.Chem import AllChem
from meeko import MoleculePreparation, PDBQTWriterLegacy
from vina import Vina

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

dock_dir = os.path.join(project_root, "01_data", "docking")
res_dir  = os.path.join(project_root, "02_results")
prep_dir = os.path.join(dock_dir, "prepared")
pose_dir = os.path.join(dock_dir, "poses")
os.makedirs(os.path.join(prep_dir, "receptors"), exist_ok=True)
os.makedirs(os.path.join(prep_dir, "ligands"), exist_ok=True)
os.makedirs(pose_dir, exist_ok=True)

SOLVENT_IONS = {
    "HOH", "DOD", "SO4", "PO4", "CL", "NA", "MG", "ZN", "CA", "MN", "FE", "K",
    "IOD", "BR", "GOL", "EDO", "DMS", "PEG", "PG4", "ACT", "TRS", "MES", "EPE",
    "MPD", "FMT", "BME", "DTT", "ACN", "ACE", "NHE", "CIT", "TAR", "MAL"
}

# 1. Prepare all Ligands via RDKit + Meeko
print("=== 1. Preparing Microbial Metabolite 3D Conformers (RDKit + Meeko) ===", flush=True)
lig_paths = {}
with open(os.path.join(res_dir, "Docking_Ligand_Manifest.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        name = row["metabolite"]
        smi  = row["smiles"]
        out_pdbqt = os.path.join(prep_dir, "ligands", f"{name}.pdbqt")
        if not os.path.exists(out_pdbqt) or os.path.getsize(out_pdbqt) < 50:
            mol = Chem.MolFromSmiles(smi)
            mol = Chem.AddHs(mol)
            params = AllChem.ETKDGv3()
            params.randomSeed = 2026
            rc = AllChem.EmbedMolecule(mol, params)
            if rc != 0:
                AllChem.EmbedMolecule(mol, useRandomCoords=True, randomSeed=2026)
            try:
                AllChem.MMFFOptimizeMolecule(mol, maxIters=500)
            except Exception:
                AllChem.UFFOptimizeMolecule(mol, maxIters=500)
            prep = MoleculePreparation()
            mol_setups = prep.prepare(mol)
            pdbqt_str, is_ok, err = PDBQTWriterLegacy.write_string(mol_setups[0])
            if not is_ok:
                raise RuntimeError(f"Meeko ligand failed for {name}: {err}")
            with open(out_pdbqt, "w", encoding="utf-8") as wf:
                wf.write(pdbqt_str)
        lig_paths[name] = out_pdbqt
        print(f"  Prepared ligand: {name}", flush=True)

# 2. Prepare Receptors + Automated Binding Box (Co-crystal HETATM pocket OR Folded Core Domain)
print("\n=== 2. Preparing Receptors & Pocket Boxes ===", flush=True)
mk_rec = (
    os.environ.get("MEEKO_RECEPTOR")
    or shutil.which("mk_prepare_receptor.py")
    or shutil.which("mk_prepare_receptor")
)
if not mk_rec:
    raise RuntimeError("Install the Meeko receptor CLI or set MEEKO_RECEPTOR.")
rec_cli = [sys.executable, mk_rec] if mk_rec.lower().endswith(".py") else [mk_rec]

def clean_and_box(gene, raw_pdb, source_str):
    lines_atom = []
    ca_coords = []
    het_by_res = {}
    chosen_chain = None
    with open(raw_pdb, "r", encoding="utf-8", errors="replace") as f:
        all_lines = f.readlines()

    # Determine primary protein chain (prefer 'A', else first ATOM chain)
    chains = [ln[21] for ln in all_lines if ln.startswith("ATOM  ") and len(ln) >= 54]
    if chains:
        chosen_chain = "A" if "A" in chains else chains[0]

    for ln in all_lines:
        if ln.startswith("ATOM  ") and len(ln) >= 66 and ln[21] == chosen_chain:
            alt = ln[16]
            if alt not in (" ", "A"):
                continue
            ln_clean = ln[:16] + " " + ln[17:]
            # Keep high-confidence folded core if AlphaFold and protein is very elongated (>65 A)
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
            bfac = float(ln[60:66])
            lines_atom.append((ln_clean, x, y, z, bfac, ln[12:16].strip()))
        elif ln.startswith("HETATM") and len(ln) >= 54 and (chosen_chain is None or ln[21] == chosen_chain):
            resn = ln[17:20].strip()
            if resn in SOLVENT_IONS:
                continue
            key = (ln[21], resn, ln[22:26].strip())
            x, y, z = float(ln[30:38]), float(ln[38:46]), float(ln[46:54])
            het_by_res.setdefault(key, []).append([x, y, z])

    # For AlphaFold full-length models with disordered tails, restrict box center/size to pLDDT>=70 structured core if available
    all_ca = np.array([[x, y, z] for _, x, y, z, _, an in lines_atom if an == "CA"])
    hi_ca  = np.array([[x, y, z] for _, x, y, z, bf, an in lines_atom if an == "CA" and bf >= 70.0])
    if "AlphaFold" in source_str and len(hi_ca) >= 80:
        kept_lines = [l for l, _, _, _, bf, _ in lines_atom if bf >= 60.0]
        core_pts = hi_ca
    else:
        kept_lines = [l for l, *_ in lines_atom]
        core_pts = all_ca

    clean_pdb = os.path.join(prep_dir, "receptors", f"{gene}_clean.pdb")
    with open(clean_pdb, "w", encoding="utf-8") as wf:
        wf.writelines(kept_lines)
        wf.write("END\n")

    # Select largest non-solvent co-crystal HETATM ligand (>= 6 heavy atoms) if present
    best_lig = None
    best_cnt = 5
    for k, pts in het_by_res.items():
        if len(pts) > best_cnt:
            best_cnt = len(pts)
            best_lig = (k, np.array(pts))

    if best_lig is not None:
        ctr = best_lig[1].mean(axis=0)
        box = [26.0, 26.0, 26.0]
        mode = f"CoCrystal_Pocket({best_lig[0][1]}, {best_cnt} atoms)"
    else:
        ctr = core_pts.mean(axis=0)
        span = core_pts.max(axis=0) - core_pts.min(axis=0)
        box = [float(min(60.0, max(28.0, s + 6.0))) for s in span]
        mode = f"FoldedCore_Box({box[0]:.0f}x{box[1]:.0f}x{box[2]:.0f})"

    out_base = os.path.join(prep_dir, "receptors", gene)
    out_pdbqt = out_base + ".pdbqt"
    if not os.path.exists(out_pdbqt):
        cmd = [*rec_cli, "--read_pdb", clean_pdb, "-o", out_base, "-p", "--default_altloc", "A", "--allow_bad_res"]
        cp = subprocess.run(cmd, capture_output=True, text=True)
        if cp.returncode != 0 or not os.path.exists(out_pdbqt):
            # Fallback with -b (blunt ends) or full AlphaFold file
            cmd2 = [*rec_cli, "--read_pdb", raw_pdb, "-o", out_base, "-p", "--default_altloc", "A", "--allow_bad_res", "-b"]
            cp2 = subprocess.run(cmd2, capture_output=True, text=True)
            if cp2.returncode != 0 or not os.path.exists(out_pdbqt):
                raise RuntimeError(f"Receptor prep failed for {gene}: {cp.stderr}\n{cp2.stderr}")

    return out_pdbqt, [float(ctr[0]), float(ctr[1]), float(ctr[2])], box, mode

rec_info = {}
with open(os.path.join(res_dir, "Docking_Protein_Manifest.csv"), "r", encoding="utf-8") as f:
    for row in csv.DictReader(f):
        g = row["gene"]
        pdb_p = os.path.join(dock_dir, "proteins", PureWindowsPath(row["pdb_path"]).name)
        try:
            pdbqt, ctr, box, mode = clean_and_box(g, pdb_p, row["source"])
            rec_info[g] = (pdbqt, ctr, box, mode, row["source"])
            print(f"  Receptor {g:8s} READY | {mode:36s} | center=({ctr[0]:.1f},{ctr[1]:.1f},{ctr[2]:.1f})", flush=True)
        except Exception as e:
            print(f"  WARN Receptor {g} skipped: {e}", flush=True)

# 3. Run AutoDock Vina Docking Matrix
# Known gutMGene2 causal/functional links for annotation
known_links = {
    ("IDO1", "Indole-3-carboxylic_acid"): "gutMGene2: inhibition (PMID:37770127)",
    ("PEMT", "Urolithin_A"): "gutMGene2: activation (PMID:30626868) + TNBC MR FDR=0.048",
    ("DAPK3", "Urolithin_A"): "gutMGene2: activation (PMID:30626868)",
    ("TGFBR1", "Trimethylamine_N-oxide"): "gutMGene2: activation (PMID:31299216) + ER- MR",
    ("CASP1", "Caffeic_acid"): "gutMGene2: inhibition (PMID:33096252)",
    ("NLRP3", "Caffeic_acid"): "gutMGene2: inhibition (PMID:33096252) + TNBC MR",
    ("FFAR2", "Butyrate"): "gutMGene2: activation (SCFA receptor) + PureTNBC MR",
    ("FFAR2", "Propionate"): "gutMGene2: activation (SCFA receptor) + PureTNBC MR",
}

print("\n=== 3. Running AutoDock Vina 1.2.7 (exhaustiveness=16, cpu=8) ===", flush=True)
results = []
for g, (rec_pq, ctr, box, mode, src) in rec_info.items():
    v = Vina(sf_name="vina", cpu=8, seed=2026, verbosity=0)
    v.set_receptor(rec_pq)
    v.compute_vina_maps(center=ctr, box_size=box)
    for lig_name, lig_pq in lig_paths.items():
        v.set_ligand_from_file(lig_pq)
        v.dock(exhaustiveness=16, n_poses=5)
        energies = v.energies(n_poses=5)
        best_aff = float(energies[0][0])
        rmsd_lb2 = float(energies[1][1]) if len(energies) > 1 else np.nan
        pose_out = os.path.join(pose_dir, f"{g}__{lig_name}_out.pdbqt")
        v.write_poses(pose_out, n_poses=3, overwrite=True)
        ann = known_links.get((g, lig_name), "")
        results.append({
            "gene": g, "source": src, "box_mode": mode,
            "metabolite": lig_name,
            "affinity_kcal_mol": best_aff,
            "second_pose_affinity": float(energies[1][0]) if len(energies) > 1 else np.nan,
            "gutMGene_evidence": ann,
            "pose_pdbqt": os.path.relpath(pose_out, project_root).replace(os.sep, "/")
        })
        if ann or best_aff <= -6.5:
            print(f"  [Vina] {g:8s} x {lig_name:34s} -> {best_aff:6.2f} kcal/mol  {ann}", flush=True)

res_csv = os.path.join(res_dir, "Docking_Vina_Affinities_Matrix.csv")
with open(res_csv, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=list(results[0].keys()))
    w.writeheader()
    w.writerows(results)
print(f"\nSAVED {len(results)} docking scores to {res_csv}", flush=True)
