"""Audit receptor origins, HETATM stripping, search-box definitions, AlphaFold pLDDT,
UniProt topology/functional sites, contact residues, and crystal-vs-AlphaFold alignments.

Outputs:
  - 02_results/JTM_upgrade/docking/receptor_plddt_and_topology.csv
  - 01_data/JTM_upgrade/docking/qc_work/part1_audit_metrics.json
"""
import csv
import json
import os
import shutil
import numpy as np
import pandas as pd
from Bio.PDB import PDBParser, MMCIFParser, Superimposer

BASE = "/mnt/d/Bioinformatics/TNBC_Microbiome"
DOCK = f"{BASE}/01_data/docking"
JTM_DATA = f"{BASE}/01_data/JTM_upgrade/docking"
RES = f"{BASE}/02_results"
OUT = f"{BASE}/02_results/JTM_upgrade/docking"
WORK = f"{JTM_DATA}/qc_work"
os.makedirs(OUT, exist_ok=True)
os.makedirs(WORK, exist_ok=True)

SOLVENT_IONS = {
    "HOH", "DOD", "SO4", "PO4", "CL", "NA", "MG", "ZN", "CA", "MN", "FE", "K",
    "IOD", "BR", "GOL", "EDO", "DMS", "PEG", "PG4", "ACT", "TRS", "MES", "EPE",
    "MPD", "FMT", "BME", "DTT", "ACN", "ACE", "NHE", "CIT", "TAR", "MAL"
}

def load_pdbqt_heavy_models(pq_path):
    models = []
    cur = []
    aff = None
    with open(pq_path, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("MODEL"):
                cur = []
                aff = None
            elif ln.startswith("REMARK VINA RESULT:"):
                aff = float(ln.split()[3])
            elif ln.startswith(("ATOM  ", "HETATM")):
                atype = ln[77:79].strip()
                if not atype.startswith("H"):
                    cur.append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])
            elif ln.startswith("ENDMDL"):
                models.append((aff, np.array(cur)))
    return models

def plddt_band(val):
    if val >= 90.0:
        return "very_high_>=90"
    if val >= 70.0:
        return "confident_70_90"
    if val >= 50.0:
        return "low_50_70"
    return "very_low_<50"

# ---------------------------------------------------------------------------
# 1. Audit all 11 targets from 21a/21b
# ---------------------------------------------------------------------------
manifest = pd.read_csv(f"{RES}/Docking_Protein_Manifest.csv")
aff_df = pd.read_csv(f"{RES}/Docking_Vina_Affinities_Matrix.csv")

targets_audit = []
for _, row in manifest.iterrows():
    gene = row["gene"]
    raw_p = str(row["pdb_path"]).replace("\\", "/").replace("D:/", "/mnt/d/")
    clean_p = f"{DOCK}/prepared/receptors/{gene}_clean.pdb"
    pq_p = f"{DOCK}/prepared/receptors/{gene}.pdbqt"

    raw_chains = set()
    raw_res = set()
    het_counts = {}
    ca_pts = []
    hi_ca_pts = []
    het_by_res = {}
    chosen_chain = None
    with open(raw_p, "r", encoding="utf-8", errors="replace") as f:
        lines = f.readlines()
    chains_seq = [ln[21] for ln in lines if ln.startswith("ATOM  ") and len(ln) >= 54]
    if chains_seq:
        chosen_chain = "A" if "A" in chains_seq else chains_seq[0]
    for ln in lines:
        if ln.startswith("ATOM  ") and len(ln) >= 66:
            raw_chains.add(ln[21])
            raw_res.add((ln[21], ln[22:26].strip(), ln[17:20].strip()))
            if ln[21] == chosen_chain and ln[16] in (" ", "A") and ln[12:16].strip() == "CA":
                pt = [float(ln[30:38]), float(ln[38:46]), float(ln[46:54])]
                bf = float(ln[60:66])
                ca_pts.append(pt)
                if bf >= 70.0:
                    hi_ca_pts.append(pt)
        elif ln.startswith("HETATM") and len(ln) >= 54:
            rn = ln[17:20].strip()
            het_counts[rn] = het_counts.get(rn, 0) + 1
            if (chosen_chain is None or ln[21] == chosen_chain) and rn not in SOLVENT_IONS:
                key = (ln[21], rn, ln[22:26].strip())
                het_by_res.setdefault(key, []).append([float(ln[30:38]), float(ln[38:46]), float(ln[46:54])])

    clean_res = set()
    with open(clean_p, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("ATOM  "):
                clean_res.add((ln[21], ln[22:26].strip(), ln[17:20].strip()))
    pq_res = set()
    pq_het = 0
    with open(pq_p, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("ATOM  "):
                pq_res.add((ln[21], ln[22:26].strip(), ln[17:20].strip()))
            elif ln.startswith("HETATM"):
                pq_het += 1

    # Reconstruct 21b box
    core_pts = np.array(hi_ca_pts) if ("AlphaFold" in str(row["source"]) and len(hi_ca_pts) >= 80) else np.array(ca_pts)
    best_lig = None
    best_cnt = 5
    for k, pts in het_by_res.items():
        if len(pts) > best_cnt:
            best_cnt = len(pts)
            best_lig = (k, np.array(pts))
    if best_lig is not None:
        ctr = best_lig[1].mean(axis=0)
        box = [26.0, 26.0, 26.0]
        box_type = f"CoCrystal_Pocket({best_lig[0][1]})"
    else:
        ctr = core_pts.mean(axis=0)
        span = core_pts.max(axis=0) - core_pts.min(axis=0)
        box = [float(min(60.0, max(28.0, s + 6.0))) for s in span]
        box_type = "FoldedCore_Box"

    sub_aff = aff_df[aff_df["gene"] == gene]
    box_mode_csv = sub_aff["box_mode"].iloc[0] if len(sub_aff) else ""
    targets_audit.append({
        "gene": gene,
        "uniprot": row["uniprot"],
        "best_pdb": "" if pd.isna(row["best_pdb"]) else str(row["best_pdb"]),
        "best_pdb_res": None if pd.isna(row["resolution"]) else float(row["resolution"]),
        "source_used": row["source"],
        "raw_chains": sorted(raw_chains),
        "chosen_chain": chosen_chain,
        "raw_res_count": len(raw_res),
        "clean_res_count": len(clean_res),
        "pdbqt_res_count": len(pq_res),
        "pdbqt_hetatm_count": pq_het,
        "raw_hetatms_non_hoh": {k: v for k, v in het_counts.items() if k != "HOH"},
        "raw_hoh_count": het_counts.get("HOH", 0),
        "box_mode_csv": box_mode_csv,
        "box_type": box_type,
        "box_center": [round(float(x), 4) for x in ctr],
        "box_size": [round(float(x), 4) for x in box],
        "box_volume_A3": round(float(box[0] * box[1] * box[2]), 1),
    })

# Check 9FOH entry JSON pdb_format_compatible
with open(f"{JTM_DATA}/9FOH_entry.json", "r", encoding="utf-8") as f:
    e9 = json.load(f)
foh_pdb_compat = e9.get("pdbx_database_status", {}).get("pdb_format_compatible")

# ---------------------------------------------------------------------------
# 2. Structural alignments of 2D0T, 4PK5, 9FOH onto IDO1_P14902_AF
# ---------------------------------------------------------------------------
pdb_parser = PDBParser(QUIET=True)
cif_parser = MMCIFParser(QUIET=True)

ido1_af_chain = pdb_parser.get_structure("IDO1_AF", f"{DOCK}/proteins/IDO1_P14902_AF.pdb")[0]["A"]
pemt_af_chain = pdb_parser.get_structure("PEMT_AF", f"{DOCK}/proteins/PEMT_Q9UBM1_AF.pdb")[0]["A"]
s2d0t_chain = pdb_parser.get_structure("2D0T", f"{JTM_DATA}/2D0T.pdb")[0]["A"]
s4pk5_chain = pdb_parser.get_structure("4PK5", f"{JTM_DATA}/4PK5.pdb")[0]["A"]
s9foh_chain = cif_parser.get_structure("9FOH", f"{JTM_DATA}/9FOH.cif")[0]["A"]

def align_chain_to_af(ref_chain, mob_chain):
    ref_map = {r.id[1]: r for r in ref_chain if r.id[0] == " " and "CA" in r}
    ref_ca, mob_ca = [], []
    for r in mob_chain:
        if r.id[0] == " " and "CA" in r and r.id[1] in ref_map:
            ref_ca.append(ref_map[r.id[1]]["CA"])
            mob_ca.append(r["CA"])
    sup = Superimposer()
    sup.set_atoms(ref_ca, mob_ca)
    sup.apply(list(mob_chain.get_atoms()))
    return len(ref_ca), float(sup.rms)

n_2d0t, rmsd_2d0t = align_chain_to_af(ido1_af_chain, s2d0t_chain)
n_4pk5, rmsd_4pk5 = align_chain_to_af(ido1_af_chain, s4pk5_chain)
n_9foh, rmsd_9foh = align_chain_to_af(ido1_af_chain, s9foh_chain)

def get_het_pts(chain, resname):
    for r in chain:
        if r.id[0] != " " and r.get_resname().strip() == resname:
            return np.array([a.get_coord() for a in r.get_atoms() if a.element != "H"])
    return None

hem_2d0t_pts = get_het_pts(s2d0t_chain, "HEM")
pim_2d0t_pts = get_het_pts(s2d0t_chain, "PIM")
hem_4pk5_pts = get_het_pts(s4pk5_chain, "HEM")
pkj_4pk5_pts = get_het_pts(s4pk5_chain, "PKJ")
a1iee_9foh_pts = get_het_pts(s9foh_chain, "A1IEE")

hem_2d0t_fe = None
for r in s2d0t_chain:
    if r.id[0] != " " and r.get_resname().strip() == "HEM" and "FE" in r:
        hem_2d0t_fe = r["FE"].get_coord()

# ---------------------------------------------------------------------------
# 3. Load docked poses and contact CSV
# ---------------------------------------------------------------------------
contacts_df = pd.read_csv(f"{RES}/Docking_Key_Complexes_Contact_Residues.csv")
pemt_pose_models = load_pdbqt_heavy_models(f"{DOCK}/poses/PEMT__Urolithin_A_out.pdbqt")
ido1_pose_models = load_pdbqt_heavy_models(f"{DOCK}/poses/IDO1__Indole-3-carboxylic_acid_out.pdbqt")

# ---------------------------------------------------------------------------
# 4. Build per-residue table for PEMT and IDO1
# ---------------------------------------------------------------------------
def parse_uniprot_annotations(json_path, seq_len, gene):
    with open(json_path, "r", encoding="utf-8") as f:
        u = json.load(f)
    topo = {i: ("Cytoplasmic" if gene == "IDO1" else "Unannotated") for i in range(1, seq_len + 1)}
    is_tm = {i: False for i in range(1, seq_len + 1)}
    reg = {i: [] for i in range(1, seq_len + 1)}
    func = {i: [] for i in range(1, seq_len + 1)}
    for feat in u.get("features", []):
        ft = feat.get("type")
        loc = feat.get("location", {})
        st = loc.get("start", {}).get("value")
        en = loc.get("end", {}).get("value")
        if st is None or en is None:
            continue
        desc = feat.get("description", "")
        lig = feat.get("ligand", {}).get("name", "") if feat.get("ligand") else ""
        if ft in ("Topological domain", "Transmembrane", "Intramembrane"):
            lbl = f"{ft}: {desc}" if desc else ft
            for i in range(int(st), int(en) + 1):
                topo[i] = lbl
                if ft in ("Transmembrane", "Intramembrane"):
                    is_tm[i] = True
        elif ft in ("Region", "Disordered", "Domain"):
            lbl = f"{ft}:{desc}" if desc else ft
            for i in range(int(st), int(en) + 1):
                if lbl not in reg[i]:
                    reg[i].append(lbl)
        elif ft in ("Binding site", "Active site", "Metal binding", "Site"):
            lbl = f"{ft} ({lig}; {desc})".strip("; )") + (")" if "(" in f"{ft} ({lig}" and not f"{ft} ({lig}; {desc})".endswith(")") else "")
            lbl = f"{ft}: {lig} [{desc}]".strip(" []") if lig else f"{ft}: {desc}"
            for i in range(int(st), int(en) + 1):
                func[i].append(lbl)
        elif ft == "Mutagenesis":
            for i in range(int(st), int(en) + 1):
                func[i].append(f"Mutagenesis: {desc}")
    return u, topo, is_tm, reg, func

def get_kept_resnums(path):
    s = set()
    with open(path, "r", encoding="utf-8") as f:
        for ln in f:
            if ln.startswith("ATOM  "):
                s.add(int(ln[22:26]))
    return s

def min_dist(pts_a, pts_b):
    if pts_a is None or pts_b is None or len(pts_a) == 0 or len(pts_b) == 0:
        return np.nan
    return float(np.linalg.norm(pts_a[:, None, :] - pts_b[None, :, :], axis=-1).min())

rows = []
summary_stats = {}

for gene, acc, af_chain, pose_models, lig_name in [
    ("PEMT", "Q9UBM1", pemt_af_chain, pemt_pose_models, "Urolithin_A"),
    ("IDO1", "P14902", ido1_af_chain, ido1_pose_models, "Indole-3-carboxylic_acid"),
]:
    res_list = [r for r in af_chain if r.id[0] == " "]
    seq_len = len(res_list)
    u, topo, is_tm, reg, func = parse_uniprot_annotations(f"{JTM_DATA}/{gene}_uniprot.json", seq_len, gene)
    kept_clean = get_kept_resnums(f"{DOCK}/prepared/receptors/{gene}_clean.pdb")
    kept_pq = get_kept_resnums(f"{DOCK}/prepared/receptors/{gene}.pdbqt")

    c_sub = contacts_df[(contacts_df["gene"] == gene) & (contacts_df["metabolite"] == lig_name)]
    c_map = {int(r["resi"]): float(r["min_dist"]) for _, r in c_sub.iterrows()}

    plddt_vals = []
    contact_plddts = []
    tm_plddts = []

    for r in res_list:
        rnum = int(r.id[1])
        rname = r.get_resname().strip()
        plddt = float(r["CA"].get_bfactor())
        band = plddt_band(plddt)
        plddt_vals.append(plddt)
        if is_tm[rnum]:
            tm_plddts.append(plddt)

        rpts = np.array([a.get_coord() for a in r.get_atoms() if a.element != "H"])
        d_p1 = min_dist(rpts, pose_models[0][1]) if len(pose_models) >= 1 else np.nan
        d_p2 = min_dist(rpts, pose_models[1][1]) if len(pose_models) >= 2 else np.nan
        d_p3 = min_dist(rpts, pose_models[2][1]) if len(pose_models) >= 3 else np.nan

        in_csv = rnum in c_map
        if in_csv:
            contact_plddts.append(plddt)

        d_hem = min_dist(rpts, hem_2d0t_pts) if gene == "IDO1" else np.nan
        d_pim = min_dist(rpts, pim_2d0t_pts) if gene == "IDO1" else np.nan
        d_a1iee = min_dist(rpts, a1iee_9foh_pts) if gene == "IDO1" else np.nan

        rows.append({
            "gene": gene,
            "uniprot_acc": acc,
            "resnum": rnum,
            "resname": rname,
            "plddt": round(plddt, 2),
            "plddt_band": band,
            "kept_in_21b_clean_pdb": rnum in kept_clean,
            "kept_in_21b_pdbqt": rnum in kept_pq,
            "uniprot_topology": topo[rnum],
            "is_transmembrane_or_intramembrane": is_tm[rnum],
            "uniprot_region_or_disordered": "; ".join(reg[rnum]),
            "uniprot_functional_site": "; ".join(func[rnum]),
            "in_docking_contacts_csv": in_csv,
            "contacts_csv_min_dist_A": round(c_map[rnum], 4) if in_csv else "",
            "pose1_min_heavy_dist_A": round(d_p1, 4),
            "pose2_min_heavy_dist_A": round(d_p2, 4),
            "pose3_min_heavy_dist_A": round(d_p3, 4),
            "dist_to_aligned_2D0T_HEM_A": round(d_hem, 4) if not np.isnan(d_hem) else "",
            "dist_to_aligned_2D0T_PIM_A": round(d_pim, 4) if not np.isnan(d_pim) else "",
            "dist_to_aligned_9FOH_A1IEE_A": round(d_a1iee, 4) if not np.isnan(d_a1iee) else "",
        })

    arr = np.array(plddt_vals)
    c_arr = np.array(contact_plddts)
    # Check pose-to-pose centroid separation
    pose1_pts = pose_models[0][1]
    pose2_pts = pose_models[1][1]
    pose_info = {
        "n_res": seq_len,
        "plddt_mean": round(float(arr.mean()), 2),
        "plddt_sd": round(float(arr.std()), 2),
        "plddt_median": round(float(np.median(arr)), 2),
        "plddt_iqr": [round(float(np.percentile(arr, 25)), 2), round(float(np.percentile(arr, 75)), 2)],
        "plddt_min": round(float(arr.min()), 2),
        "plddt_max": round(float(arr.max()), 2),
        "n_very_high_ge90": int((arr >= 90).sum()),
        "n_confident_70_90": int(((arr >= 70) & (arr < 90)).sum()),
        "n_low_50_70": int(((arr >= 50) & (arr < 70)).sum()),
        "n_very_low_lt50": int((arr < 50).sum()),
        "n_dropped_by_21b_lt60": int(seq_len - len(kept_clean)),
        "dropped_resnums": sorted(set(range(1, seq_len + 1)) - kept_clean),
        "n_contacts_csv": len(c_map),
        "contacts_csv_resnums": sorted(c_map.keys()),
        "contacts_plddt_mean": round(float(c_arr.mean()), 2),
        "contacts_plddt_min": round(float(c_arr.min()), 2),
        "contacts_plddt_max": round(float(c_arr.max()), 2),
        "pose1_aff": pose_models[0][0],
        "pose2_aff": pose_models[1][0],
        "pose1_2_delta_aff": round(float(pose_models[1][0] - pose_models[0][0]), 3),
        "pose1_2_centroid_shift_A": round(float(np.linalg.norm(pose1_pts.mean(0) - pose2_pts.mean(0))), 3),
    }
    if tm_plddts:
        tm_arr = np.array(tm_plddts)
        pose_info["tm_res_count"] = len(tm_arr)
        pose_info["tm_plddt_mean"] = round(float(tm_arr.mean()), 2)
        pose_info["tm_plddt_min"] = round(float(tm_arr.min()), 2)
    summary_stats[gene] = pose_info

# Add IDO1 alignment overlap metrics
ido1_p1 = ido1_pose_models[0][1]
summary_stats["IDO1_crystal_overlap"] = {
    "2D0T_aligned_CA_count": n_2d0t,
    "2D0T_CA_RMSD_A": round(rmsd_2d0t, 3),
    "2D0T_HEM_min_dist_to_I3CA_pose1_A": round(min_dist(ido1_p1, hem_2d0t_pts), 3),
    "2D0T_HEM_FE_dist_to_I3CA_pose1_min_A": round(float(np.linalg.norm(ido1_p1 - hem_2d0t_fe, axis=1).min()), 3),
    "2D0T_PIM_min_dist_to_I3CA_pose1_A": round(min_dist(ido1_p1, pim_2d0t_pts), 3),
    "4PK5_aligned_CA_count": n_4pk5,
    "4PK5_CA_RMSD_A": round(rmsd_4pk5, 3),
    "4PK5_HEM_min_dist_to_I3CA_pose1_A": round(min_dist(ido1_p1, hem_4pk5_pts), 3),
    "4PK5_PKJ_min_dist_to_I3CA_pose1_A": round(min_dist(ido1_p1, pkj_4pk5_pts), 3),
    "9FOH_aligned_CA_count": n_9foh,
    "9FOH_CA_RMSD_A": round(rmsd_9foh, 3),
    "9FOH_A1IEE_min_dist_to_I3CA_pose1_A": round(min_dist(ido1_p1, a1iee_9foh_pts), 3),
    "9FOH_pdb_format_compatible": foh_pdb_compat,
    "HIS346_NE2_to_I3CA_pose1_min_A": round(float(np.linalg.norm(ido1_p1 - ido1_af_chain[346]["NE2"].get_coord(), axis=1).min()), 3),
}

out_df = pd.DataFrame(rows)
out_csv = f"{OUT}/receptor_plddt_and_topology.csv"
out_df.to_csv(out_csv, index=False)

# Also record environment tool inventory
env_info = {
    "wsl_python_exe": shutil.which("python3") or "/root/miniconda3/envs/bioinfo-singlecell-python/bin/python",
    "cli_tools": {
        cmd: shutil.which(cmd) for cmd in ["vina", "obabel", "autogrid4", "autodock4", "prepare_receptor4.py", "mk_prepare_receptor.py", "mk_prepare_ligand.py"]
    }
}
import vina, meeko, rdkit, Bio
env_info["python_packages"] = {
    "vina": getattr(vina, "__version__", "unknown"),
    "meeko": getattr(meeko, "__version__", "unknown"),
    "rdkit": rdkit.__version__,
    "biopython": Bio.__version__,
}

metrics = {
    "environment": env_info,
    "targets_audit": targets_audit,
    "summary_stats": summary_stats,
}
with open(f"{WORK}/part1_audit_metrics.json", "w", encoding="utf-8") as f:
    json.dump(metrics, f, indent=2)

print(f"Saved {out_csv} ({len(out_df)} rows)", flush=True)
print(json.dumps(metrics["summary_stats"], indent=2), flush=True)
print("9FOH pdb_format_compatible:", foh_pdb_compat, flush=True)
