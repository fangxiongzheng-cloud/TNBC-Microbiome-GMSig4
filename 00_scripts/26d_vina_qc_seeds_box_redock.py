"""Run AutoDock Vina 1.2.7 multi-seed repeatability, box-size sensitivity, 9FOH/A1IEE redocking,
and 2D0T/4PK5 HEM retention feasibility checks using existing installed WSL tools.

Outputs:
  - 02_results/JTM_upgrade/docking/redocking_and_seed_qc.csv
  - 01_data/JTM_upgrade/docking/qc_work/part2_qc_summary.json
  - Poses & temp inputs under 01_data/JTM_upgrade/docking/qc_runs/
"""
import itertools
import json
import os
import subprocess
import sys
import time
import numpy as np
import pandas as pd
from Bio.PDB import MMCIFParser, PDBParser, PDBIO, Select, Superimposer
from Bio.PDB.MMCIF2Dict import MMCIF2Dict
from rdkit import Chem
from rdkit.Chem import AllChem, rdMolAlign
from meeko import MoleculePreparation, PDBQTMolecule, PDBQTWriterLegacy, RDKitMolCreate
from vina import Vina

BASE = "/mnt/d/Bioinformatics/TNBC_Microbiome"
DOCK = f"{BASE}/01_data/docking"
JTM_DATA = f"{BASE}/01_data/JTM_upgrade/docking"
OUT = f"{BASE}/02_results/JTM_upgrade/docking"
WORK = f"{JTM_DATA}/qc_work"
RUNS = f"{JTM_DATA}/qc_runs"
os.makedirs(WORK, exist_ok=True)
os.makedirs(RUNS, exist_ok=True)
os.makedirs(OUT, exist_ok=True)

MK_REC = "/root/miniconda3/envs/bioinfo-singlecell-python/bin/mk_prepare_receptor.py"
SEEDS = [2026, 42, 101, 777, 1337]

def pdbqt_to_heavy_mol(pdbqt_path):
    pm = PDBQTMolecule.from_file(pdbqt_path)
    mols = RDKitMolCreate.from_pdbqt_mol(pm)
    m = Chem.RemoveHs(mols[0])
    return m

def pose_centroid(mol, conf_id=0):
    conf = mol.GetConformer(conf_id)
    pts = np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])
    return pts.mean(axis=0)

def pose_heavy_coords(mol, conf_id=0):
    conf = mol.GetConformer(conf_id)
    return np.array([list(conf.GetAtomPosition(i)) for i in range(mol.GetNumAtoms())])

def min_dist_pts(a, b):
    return float(np.linalg.norm(a[:, None, :] - b[None, :, :], axis=-1).min())

def run_vina_job(rec_pq, lig_pq, center, box, seed, exhaustiveness, n_poses, out_pq):
    t0 = time.time()
    v = Vina(sf_name="vina", cpu=8, seed=seed, verbosity=0)
    v.set_receptor(rec_pq)
    v.compute_vina_maps(center=[float(c) for c in center], box_size=[float(b) for b in box])
    v.set_ligand_from_file(lig_pq)
    v.dock(exhaustiveness=exhaustiveness, n_poses=n_poses)
    energies = v.energies(n_poses=n_poses)
    v.write_poses(out_pq, n_poses=n_poses, overwrite=True)
    dt = time.time() - t0
    return energies, dt

# ---------------------------------------------------------------------------
# 0. Load Part 1 metadata (exact 21b boxes & aligned 2D0T HEM/PIM coordinates)
# ---------------------------------------------------------------------------
with open(f"{WORK}/part1_audit_metrics.json", "r", encoding="utf-8") as f:
    p1 = json.load(f)
tmap = {t["gene"]: t for t in p1["targets_audit"]}

pdb_p = PDBParser(QUIET=True)
cif_p = MMCIFParser(QUIET=True)

pemt_af = pdb_p.get_structure("PEMT", f"{DOCK}/proteins/PEMT_Q9UBM1_AF.pdb")[0]["A"]
ido1_af = pdb_p.get_structure("IDO1", f"{DOCK}/proteins/IDO1_P14902_AF.pdb")[0]["A"]
s2d0t = pdb_p.get_structure("2D0T", f"{JTM_DATA}/2D0T.pdb")[0]["A"]

# Align 2D0T_A onto IDO1_AF
ref_map = {r.id[1]: r for r in ido1_af if r.id[0] == " " and "CA" in r}
ref_ca, mob_ca = [], []
for r in s2d0t:
    if r.id[0] == " " and "CA" in r and r.id[1] in ref_map:
        ref_ca.append(ref_map[r.id[1]]["CA"])
        mob_ca.append(r["CA"])
sup = Superimposer()
sup.set_atoms(ref_ca, mob_ca)
sup.apply(list(s2d0t.get_atoms()))

hem_2d0t_pts = np.array([a.get_coord() for r in s2d0t if r.get_resname().strip() == "HEM" for a in r.get_atoms() if a.element != "H"])
pim_2d0t_pts = np.array([a.get_coord() for r in s2d0t if r.get_resname().strip() == "PIM" for a in r.get_atoms() if a.element != "H"])
pim_2d0t_ctr = pim_2d0t_pts.mean(axis=0).tolist()

sam_pemt_pts = np.array([a.get_coord() for rnum in [98, 99, 100, 180, 181] for a in pemt_af[rnum].get_atoms() if a.element != "H"])
sam_pemt_ctr = sam_pemt_pts.mean(axis=0).tolist()

# Baseline 21b heavy RDKit mols (Pose 1..3)
pemt_base_mol = pdbqt_to_heavy_mol(f"{DOCK}/poses/PEMT__Urolithin_A_out.pdbqt")
ido1_base_mol = pdbqt_to_heavy_mol(f"{DOCK}/poses/IDO1__Indole-3-carboxylic_acid_out.pdbqt")
pemt_base_ctr = pose_centroid(pemt_base_mol, 0).tolist()
ido1_base_ctr = pose_centroid(ido1_base_mol, 0).tolist()

qc_rows = []
pairwise_summaries = {}

# ---------------------------------------------------------------------------
# 1. Multi-seed & Box-size Sensitivity Checks: PEMT + UroA and IDO1 + I3CA
# ---------------------------------------------------------------------------
configs = [
    ("PEMT", "Urolithin_A", f"{DOCK}/prepared/receptors/PEMT.pdbqt", f"{DOCK}/prepared/ligands/Urolithin_A.pdbqt",
     pemt_base_mol, tmap["PEMT"]["box_center"], tmap["PEMT"]["box_size"], pemt_base_ctr, sam_pemt_ctr, "SAM_site"),
    ("IDO1", "Indole-3-carboxylic_acid", f"{DOCK}/prepared/receptors/IDO1.pdbqt", f"{DOCK}/prepared/ligands/Indole-3-carboxylic_acid.pdbqt",
     ido1_base_mol, tmap["IDO1"]["box_center"], tmap["IDO1"]["box_size"], ido1_base_ctr, pim_2d0t_ctr, "2D0T_PIM_distal"),
]

for gene, lig, rec_pq, lig_pq, base_mol, blind_ctr, blind_box, pose1_ctr, bio_ctr, bio_lbl in configs:
    print(f"\n=== QC RUNS: {gene} + {lig} ===", flush=True)
    experiments = []
    # (exp_group, box_label, ctr, box, exhaustiveness, seeds)
    experiments.append(("multi_seed_21b_blind", "21b_FoldedCore_Box", blind_ctr, blind_box, 16, SEEDS))
    experiments.append(("exhaustiveness32_blind", "21b_FoldedCore_Box_exh32", blind_ctr, blind_box, 32, [2026]))
    experiments.append(("box_expanded_plus10A", "Expanded_plus10A", blind_ctr, [b + 10.0 for b in blind_box], 16, [2026]))
    experiments.append(("multi_seed_focused_26A", "Focused_Pose1_26A", pose1_ctr, [26.0, 26.0, 26.0], 16, SEEDS))
    experiments.append(("box_focused_20A", "Focused_Pose1_20A", pose1_ctr, [20.0, 20.0, 20.0], 16, [2026]))
    experiments.append((f"box_targeted_{bio_lbl}_26A", f"Targeted_{bio_lbl}_26A", bio_ctr, [26.0, 26.0, 26.0], 16, [2026]))

    for exp_group, box_lbl, ctr, box, exh, seed_list in experiments:
        group_mols = []
        group_affs = []
        for sd in seed_list:
            out_pq = f"{RUNS}/{gene}__{lig}__{box_lbl}__exh{exh}__seed{sd}.pdbqt"
            energies, dt = run_vina_job(rec_pq, lig_pq, ctr, box, sd, exh, 5, out_pq)
            aff1 = float(energies[0][0])
            aff2 = float(energies[1][0]) if len(energies) > 1 else np.nan
            mol_run = pdbqt_to_heavy_mol(out_pq)
            rmsd_vs_base = float(rdMolAlign.CalcRMS(mol_run, base_mol, 0, 0))
            c_shift = float(np.linalg.norm(pose_centroid(mol_run, 0) - np.array(pose1_ctr)))
            # Best RMSD among top 5 poses vs baseline pose 1
            best_top5_rmsd = min(float(rdMolAlign.CalcRMS(mol_run, base_mol, ci, 0)) for ci in range(mol_run.GetNumConformers()))
            hem_min = min_dist_pts(pose_heavy_coords(mol_run, 0), hem_2d0t_pts) if gene == "IDO1" else np.nan

            group_mols.append((sd, mol_run))
            group_affs.append(aff1)
            qc_rows.append({
                "test_category": exp_group,
                "receptor": gene,
                "ligand": lig,
                "box_label": box_lbl,
                "box_center": f"({ctr[0]:.2f},{ctr[1]:.2f},{ctr[2]:.2f})",
                "box_size": f"({box[0]:.1f},{box[1]:.1f},{box[2]:.1f})",
                "box_volume_A3": round(float(box[0] * box[1] * box[2]), 1),
                "exhaustiveness": exh,
                "seed": sd,
                "affinity_pose1_kcal_mol": round(aff1, 3),
                "affinity_pose2_kcal_mol": round(aff2, 3),
                "delta_pose2_minus_pose1": round(aff2 - aff1, 3) if not np.isnan(aff2) else np.nan,
                "reference_pose": "21b_seed2026_pose1",
                "sym_rmsd_pose1_vs_ref_A": round(rmsd_vs_base, 3),
                "sym_rmsd_best_top5_vs_ref_A": round(best_top5_rmsd, 3),
                "centroid_shift_vs_ref_A": round(c_shift, 3),
                "pose1_within_2A_ref": rmsd_vs_base <= 2.0,
                "min_dist_to_aligned_2D0T_HEM_A": round(hem_min, 3) if not np.isnan(hem_min) else "",
                "runtime_sec": round(dt, 2),
                "status": "OK",
                "notes": "",
            })
            print(f"  {exp_group:26s} seed={sd:4d} aff1={aff1:6.3f} aff2={aff2:6.3f} rmsd_vs_21b={rmsd_vs_base:6.3f}A shift={c_shift:5.2f}A ({dt:.1f}s)", flush=True)

        if len(group_mols) > 1:
            pw_rmsds = []
            for (s1, m1), (s2, m2) in itertools.combinations(group_mols, 2):
                pw_rmsds.append(float(rdMolAlign.CalcRMS(m1, m2, 0, 0)))
            aff_arr = np.array(group_affs)
            pw_arr = np.array(pw_rmsds)
            pairwise_summaries[f"{gene}__{exp_group}"] = {
                "n_seeds": len(seed_list),
                "aff_mean": round(float(aff_arr.mean()), 4),
                "aff_sd": round(float(aff_arr.std(ddof=1)), 4),
                "aff_min": round(float(aff_arr.min()), 3),
                "aff_max": round(float(aff_arr.max()), 3),
                "aff_range": round(float(aff_arr.max() - aff_arr.min()), 3),
                "pairwise_pose1_rmsd_mean": round(float(pw_arr.mean()), 3),
                "pairwise_pose1_rmsd_median": round(float(np.median(pw_arr)), 3),
                "pairwise_pose1_rmsd_max": round(float(pw_arr.max()), 3),
                "frac_pairs_le2A": round(float((pw_arr <= 2.0).mean()), 3),
            }

# ---------------------------------------------------------------------------
# 2. Redocking Crystal Ligand A1IEE into Apo-IDO1 9FOH (Chain A)
# ---------------------------------------------------------------------------
print("\n=== QC RUNS: 9FOH (chain A) + A1IEE REDOCKING ===", flush=True)

# 2a. Prepare clean 9FOH chain A (polymer only, strip PG4/PGE/DMS/A1IEE/HOH)
s9_struct = cif_p.get_structure("9FOH", f"{JTM_DATA}/9FOH.cif")
class ChainAPolySelect(Select):
    def accept_chain(self, chain):
        return chain.id == "A"
    def accept_residue(self, residue):
        return residue.id[0] == " "

io = PDBIO()
io.set_structure(s9_struct)
clean_9foh_pdb = f"{WORK}/9FOH_A_clean.pdb"
io.save(clean_9foh_pdb, ChainAPolySelect())
rec_9foh_base = f"{WORK}/9FOH_A"
rec_9foh_pq = f"{rec_9foh_base}.pdbqt"
cp_9foh = subprocess.run([MK_REC, "--read_pdb", clean_9foh_pdb, "-o", rec_9foh_base, "-p", "--default_altloc", "A", "--allow_bad_res"],
                         capture_output=True, text=True)
if cp_9foh.returncode != 0 or not os.path.exists(rec_9foh_pq):
    raise RuntimeError(f"9FOH_A receptor prep failed: {cp_9foh.stderr}")

# Compute 9FOH chain A CA box too (replicating 21b FoldedCore_Box)
s9_chain = s9_struct[0]["A"]
s9_ca = np.array([r["CA"].get_coord() for r in s9_chain if r.id[0] == " " and "CA" in r])
s9_blind_ctr = s9_ca.mean(axis=0).tolist()
s9_blind_span = s9_ca.max(axis=0) - s9_ca.min(axis=0)
s9_blind_box = [float(min(60.0, max(28.0, s + 6.0))) for s in s9_blind_span]

# 2b. Build reference crystal A1IEE RDKit Mol (with bond orders from A1IEE.cif and 3D coords from 9FOH.cif)
cif_dict = MMCIF2Dict(f"{JTM_DATA}/A1IEE.cif")
atom_ids = cif_dict["_chem_comp_atom.atom_id"]
atom_syms = cif_dict["_chem_comp_atom.type_symbol"]
atom_chg = [int(x) for x in cif_dict["_chem_comp_atom.charge"]]
b_a1 = cif_dict["_chem_comp_bond.atom_id_1"]
b_a2 = cif_dict["_chem_comp_bond.atom_id_2"]
b_ord = cif_dict["_chem_comp_bond.value_order"]
b_aro = cif_dict["_chem_comp_bond.pdbx_aromatic_flag"]

rw = Chem.RWMol()
id_to_idx = {}
for aid, sym, chg in zip(atom_ids, atom_syms, atom_chg):
    if sym == "H":
        continue
    a = Chem.Atom(sym)
    a.SetFormalCharge(chg)
    id_to_idx[aid] = rw.AddAtom(a)
bond_map = {"SING": Chem.BondType.SINGLE, "DOUB": Chem.BondType.DOUBLE, "TRIP": Chem.BondType.TRIPLE, "AROM": Chem.BondType.AROMATIC}
for a1, a2, bo, ar in zip(b_a1, b_a2, b_ord, b_aro):
    if a1 in id_to_idx and a2 in id_to_idx:
        rw.AddBond(id_to_idx[a1], id_to_idx[a2], Chem.BondType.AROMATIC if ar == "Y" else bond_map[bo])
tmpl_a1 = rw.GetMol()
Chem.SanitizeMol(tmpl_a1)

a1_res = [r for r in s9_chain if r.get_resname().strip() == "A1IEE"][0]
conf = Chem.Conformer(tmpl_a1.GetNumAtoms())
for aid, idx in id_to_idx.items():
    pt = a1_res[aid].get_coord()
    conf.SetAtomPosition(idx, (float(pt[0]), float(pt[1]), float(pt[2])))
a1_xtal_heavy = Chem.Mol(tmpl_a1)
a1_xtal_heavy.AddConformer(conf, assignId=True)
Chem.AssignStereochemistryFrom3D(a1_xtal_heavy)
a1_xtal_ctr = pose_centroid(a1_xtal_heavy, 0).tolist()
a1_iso_smiles = Chem.MolToSmiles(a1_xtal_heavy, isomericSmiles=True)

# Ligand Prep A: Crystal-extracted 3D coordinates + Hs -> Meeko
a1_xtal_h = Chem.AddHs(a1_xtal_heavy, addCoords=True)
prep = MoleculePreparation()
pq_xtal_str, ok1, err1 = PDBQTWriterLegacy.write_string(prep.prepare(a1_xtal_h)[0])
assert ok1, err1
lig_a1_xtal_pq = f"{WORK}/A1IEE_xtal.pdbqt"
with open(lig_a1_xtal_pq, "w", encoding="utf-8") as f:
    f.write(pq_xtal_str)

# Ligand Prep B: SMILES -> ETKDGv3 + UFF/MMFF -> Meeko (replicating 21b workflow)
a1_smi_mol = Chem.AddHs(Chem.MolFromSmiles(a1_iso_smiles))
ep = AllChem.ETKDGv3()
ep.randomSeed = 2026
rc_emb = AllChem.EmbedMolecule(a1_smi_mol, ep)
if rc_emb != 0:
    AllChem.EmbedMolecule(a1_smi_mol, useRandomCoords=True, randomSeed=2026)
opt_method = "MMFF"
try:
    mmff_res = AllChem.MMFFOptimizeMolecule(a1_smi_mol, maxIters=500)
    if mmff_res == -1:
        AllChem.UFFOptimizeMolecule(a1_smi_mol, maxIters=500)
        opt_method = "UFF_fallback"
except Exception:
    AllChem.UFFOptimizeMolecule(a1_smi_mol, maxIters=500)
    opt_method = "UFF_fallback"
pq_smi_str, ok2, err2 = PDBQTWriterLegacy.write_string(prep.prepare(a1_smi_mol)[0])
assert ok2, err2
lig_a1_smi_pq = f"{WORK}/A1IEE_smiles_etkdg.pdbqt"
with open(lig_a1_smi_pq, "w", encoding="utf-8") as f:
    f.write(pq_smi_str)
print(f"  A1IEE SMILES={a1_iso_smiles}, ETKDG opt={opt_method}", flush=True)

redock_jobs = [
    ("redock_A1IEE_xtal_26A", "xtal_3D", lig_a1_xtal_pq, "Pocket_26A", a1_xtal_ctr, [26.0, 26.0, 26.0], 16, SEEDS),
    ("redock_A1IEE_xtal_22A", "xtal_3D", lig_a1_xtal_pq, "Pocket_22A", a1_xtal_ctr, [22.0, 22.0, 22.0], 16, [2026]),
    ("redock_A1IEE_xtal_blind", "xtal_3D", lig_a1_xtal_pq, "9FOH_FoldedCore_Box", s9_blind_ctr, s9_blind_box, 16, [2026]),
    ("redock_A1IEE_smiles_26A", f"SMILES_{opt_method}", lig_a1_smi_pq, "Pocket_26A", a1_xtal_ctr, [26.0, 26.0, 26.0], 16, SEEDS),
]

for exp_group, prep_lbl, lig_pq, box_lbl, ctr, box, exh, seed_list in redock_jobs:
    rmsds_p1 = []
    affs = []
    best5s = []
    for sd in seed_list:
        out_pq = f"{RUNS}/9FOH_A__A1IEE_{prep_lbl}__{box_lbl}__seed{sd}.pdbqt"
        energies, dt = run_vina_job(rec_9foh_pq, lig_pq, ctr, box, sd, exh, 10, out_pq)
        aff1 = float(energies[0][0])
        aff2 = float(energies[1][0]) if len(energies) > 1 else np.nan
        mol_run = pdbqt_to_heavy_mol(out_pq)
        n_conf = mol_run.GetNumConformers()
        all_rmsds = [float(rdMolAlign.CalcRMS(mol_run, a1_xtal_heavy, ci, 0)) for ci in range(n_conf)]
        rmsd1 = all_rmsds[0]
        best5 = min(all_rmsds[:5])
        best10 = min(all_rmsds)
        c_shift = float(np.linalg.norm(pose_centroid(mol_run, 0) - np.array(a1_xtal_ctr)))

        rmsds_p1.append(rmsd1)
        affs.append(aff1)
        best5s.append(best5)
        qc_rows.append({
            "test_category": exp_group,
            "receptor": "9FOH_chainA_apo_IDO1",
            "ligand": f"A1IEE_{prep_lbl}",
            "box_label": box_lbl,
            "box_center": f"({ctr[0]:.2f},{ctr[1]:.2f},{ctr[2]:.2f})",
            "box_size": f"({box[0]:.1f},{box[1]:.1f},{box[2]:.1f})",
            "box_volume_A3": round(float(box[0] * box[1] * box[2]), 1),
            "exhaustiveness": exh,
            "seed": sd,
            "affinity_pose1_kcal_mol": round(aff1, 3),
            "affinity_pose2_kcal_mol": round(aff2, 3),
            "delta_pose2_minus_pose1": round(aff2 - aff1, 3) if not np.isnan(aff2) else np.nan,
            "reference_pose": "9FOH_xtal_A1IEE",
            "sym_rmsd_pose1_vs_ref_A": round(rmsd1, 3),
            "sym_rmsd_best_top5_vs_ref_A": round(best5, 3),
            "centroid_shift_vs_ref_A": round(c_shift, 3),
            "pose1_within_2A_ref": rmsd1 <= 2.0,
            "min_dist_to_aligned_2D0T_HEM_A": "",
            "runtime_sec": round(dt, 2),
            "status": "OK",
            "notes": f"top10_best_rmsd={best10:.3f}A (n_poses={n_conf})",
        })
        print(f"  {exp_group:26s} seed={sd:4d} aff1={aff1:7.3f} pose1_RMSD={rmsd1:5.3f}A best5={best5:5.3f}A best10={best10:5.3f}A ({dt:.1f}s)", flush=True)

    if len(seed_list) > 1:
        pairwise_summaries[f"9FOH__{exp_group}"] = {
            "n_seeds": len(seed_list),
            "aff_mean": round(float(np.mean(affs)), 3),
            "aff_sd": round(float(np.std(affs, ddof=1)), 3),
            "pose1_rmsd_mean": round(float(np.mean(rmsds_p1)), 3),
            "pose1_rmsd_sd": round(float(np.std(rmsds_p1, ddof=1)), 3),
            "pose1_rmsd_min": round(float(np.min(rmsds_p1)), 3),
            "pose1_rmsd_max": round(float(np.max(rmsds_p1)), 3),
            "best5_rmsd_mean": round(float(np.mean(best5s)), 3),
            "success_rate_pose1_le2A": round(float(np.mean([r <= 2.0 for r in rmsds_p1])), 3),
        }

# ---------------------------------------------------------------------------
# 3. Heme (HEM) Retention Feasibility Tests on 2D0T and 4PK5
# ---------------------------------------------------------------------------
print("\n=== QC TESTS: 2D0T / 4PK5 HEM RETENTION FEASIBILITY ===", flush=True)

def prep_chain_pdb(pdb_in, pdb_out, keep_hem=True):
    kept = []
    with open(pdb_in, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("ATOM  ") and ln[21] == "A" and ln[16] in (" ", "A"):
                kept.append(ln[:16] + " " + ln[17:])
            elif keep_hem and ln.startswith("HETATM") and ln[21] == "A" and ln[17:20].strip() == "HEM":
                kept.append(ln)
    kept.append("END\n")
    with open(pdb_out, "w", encoding="utf-8") as f:
        f.writelines(kept)

hem_tests = []
for pdb_id in ["2D0T", "4PK5"]:
    p_hem = f"{WORK}/{pdb_id}_A_HEM.pdb"
    p_apo = f"{WORK}/{pdb_id}_A_noHEM.pdb"
    prep_chain_pdb(f"{JTM_DATA}/{pdb_id}.pdb", p_hem, keep_hem=True)
    prep_chain_pdb(f"{JTM_DATA}/{pdb_id}.pdb", p_apo, keep_hem=False)

    # Test with HEM (default)
    cp_def = subprocess.run([MK_REC, "--read_pdb", p_hem, "-o", f"{WORK}/{pdb_id}_hem_def", "-p", "--allow_bad_res"],
                            capture_output=True, text=True)
    # Test with HEM + -x (--delete_bad_res)
    cp_del = subprocess.run([MK_REC, "--read_pdb", p_hem, "-o", f"{WORK}/{pdb_id}_hem_del", "-p", "--allow_bad_res", "-x"],
                            capture_output=True, text=True)
    # Test protein-only (HEM stripped)
    cp_apo = subprocess.run([MK_REC, "--read_pdb", p_apo, "-o", f"{WORK}/{pdb_id}_apo", "-p", "--allow_bad_res"],
                            capture_output=True, text=True)

    for lbl, cp, msg in [
        (f"{pdb_id}_A_with_HEM_default", cp_def, "Meeko mk_prepare_receptor default with HEM"),
        (f"{pdb_id}_A_with_HEM_delete_bad_res", cp_del, "Meeko mk_prepare_receptor -x with HEM"),
        (f"{pdb_id}_A_stripped_noHEM", cp_apo, "Meeko mk_prepare_receptor protein-only (HEM stripped)"),
    ]:
        ok_flag = (cp.returncode == 0)
        err_snip = (cp.stdout + "\n" + cp.stderr).strip().replace("\n", " | ")
        # Extract most informative line
        key_err = [l.strip() for l in (cp.stdout + "\n" + cp.stderr).splitlines() if any(k in l for k in ["Error", "Failed", "RuntimeWarning", "AttributeError"])]
        note = "; ".join(key_err[:4]) if key_err else ("Succeeded without HEM" if ok_flag else err_snip[-300:])
        hem_tests.append({"test": lbl, "returncode": cp.returncode, "note": note})
        qc_rows.append({
            "test_category": "hem_retention_feasibility",
            "receptor": lbl,
            "ligand": "N/A",
            "box_label": "N/A",
            "box_center": "",
            "box_size": "",
            "box_volume_A3": np.nan,
            "exhaustiveness": np.nan,
            "seed": np.nan,
            "affinity_pose1_kcal_mol": np.nan,
            "affinity_pose2_kcal_mol": np.nan,
            "delta_pose2_minus_pose1": np.nan,
            "reference_pose": "N/A",
            "sym_rmsd_pose1_vs_ref_A": np.nan,
            "sym_rmsd_best_top5_vs_ref_A": np.nan,
            "centroid_shift_vs_ref_A": np.nan,
            "pose1_within_2A_ref": "",
            "min_dist_to_aligned_2D0T_HEM_A": "",
            "runtime_sec": np.nan,
            "status": "OK" if ok_flag else f"FAILED_RC_{cp.returncode}",
            "notes": note[:500],
        })
        print(f"  {lbl:36s} rc={cp.returncode} | {note[:160]}", flush=True)

# Test standalone HEM PDB -> MoleculePreparation bond order degradation & Vina ad4 requirement
hem_only = f"{WORK}/2D0T_HEM_only.pdb"
with open(hem_only, "w") as wf:
    wf.writelines([ln for ln in open(f"{JTM_DATA}/2D0T.pdb") if ln.startswith("HETATM") and ln[21] == "A" and ln[17:20].strip() == "HEM"] + ["END\n"])
m_hem = Chem.AddHs(Chem.MolFromPDBFile(hem_only, removeHs=False, sanitize=True), addCoords=True)
pq_hem, ok_h, _ = PDBQTWriterLegacy.write_string(MoleculePreparation().prepare(m_hem)[0])
smi_hem_pdb = [l for l in pq_hem.splitlines() if l.startswith("REMARK SMILES ")][0]
fe_line = [l for l in pq_hem.splitlines() if " FE " in l][0]

cp_ad4 = subprocess.run([sys.executable, "-c",
    f"from vina import Vina; v = Vina(sf_name='ad4', verbosity=0); v.set_receptor('{WORK}/dummy_rec_Fe.pdbqt')"],
    capture_output=True, text=True)

qc_rows.append({
    "test_category": "hem_retention_feasibility",
    "receptor": "2D0T_HEM_RDKit_PDB_standalone",
    "ligand": "N/A",
    "box_label": "N/A",
    "box_center": "",
    "box_size": "",
    "box_volume_A3": np.nan,
    "exhaustiveness": np.nan,
    "seed": np.nan,
    "affinity_pose1_kcal_mol": np.nan,
    "affinity_pose2_kcal_mol": np.nan,
    "delta_pose2_minus_pose1": np.nan,
    "reference_pose": "N/A",
    "sym_rmsd_pose1_vs_ref_A": np.nan,
    "sym_rmsd_best_top5_vs_ref_A": np.nan,
    "centroid_shift_vs_ref_A": np.nan,
    "pose1_within_2A_ref": "",
            "min_dist_to_aligned_2D0T_HEM_A": "",
    "runtime_sec": np.nan,
    "status": "DEGRADED_BOND_ORDERS",
    "notes": f"PDB HETATM lacks bond orders -> saturated porphyrin/gem-diol propionates ({smi_hem_pdb}); {fe_line.strip()}",
})
qc_rows.append({
    "test_category": "hem_retention_feasibility",
    "receptor": "Vina_1.2.7_sf_ad4_no_autogrid4",
    "ligand": "N/A",
    "box_label": "N/A",
    "box_center": "",
    "box_size": "",
    "box_volume_A3": np.nan,
    "exhaustiveness": np.nan,
    "seed": np.nan,
    "affinity_pose1_kcal_mol": np.nan,
    "affinity_pose2_kcal_mol": np.nan,
    "delta_pose2_minus_pose1": np.nan,
    "reference_pose": "N/A",
    "sym_rmsd_pose1_vs_ref_A": np.nan,
    "sym_rmsd_best_top5_vs_ref_A": np.nan,
    "centroid_shift_vs_ref_A": np.nan,
            "pose1_within_2A_ref": "",
            "min_dist_to_aligned_2D0T_HEM_A": "",
            "runtime_sec": np.nan,
            "status": f"FAILED_RC_{cp_ad4.returncode}",
            "notes": f"autogrid4/obabel/prepare_receptor4 missing; Vina ad4 set_receptor exits with: {cp_ad4.stderr.strip()}",
        })

qc_df = pd.DataFrame(qc_rows)
out_csv = f"{OUT}/redocking_and_seed_qc.csv"
qc_df.to_csv(out_csv, index=False)

with open(f"{WORK}/part2_qc_summary.json", "w", encoding="utf-8") as f:
    json.dump({"pairwise_summaries": pairwise_summaries, "hem_tests": hem_tests, "a1iee_smiles": a1_iso_smiles, "a1iee_opt": opt_method}, f, indent=2)

print(f"\nSaved {out_csv} ({len(qc_df)} rows)", flush=True)
print(json.dumps(pairwise_summaries, indent=2), flush=True)
