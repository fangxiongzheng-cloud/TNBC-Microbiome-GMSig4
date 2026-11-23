"""Stage C1: Pre-specified non-circular patient-level spatial transcriptomics analysis
on GSE210616 (Bassiouni et al., Cancer Res 2023; N=22 primary TNBC patients, 43 Visium sections).
Uses exact frozen GMSig-4 Ridge weights and GMSig-4-free compartment/CD8 gene sets."""

import gzip
import json
from pathlib import Path
import h5py
import numpy as np
import pandas as pd
from scipy.sparse import csc_matrix
from scipy.stats import spearmanr, mannwhitneyu

ROOT = Path("/mnt/d/Bioinformatics/TNBC_Microbiome")
DATA = ROOT / "01_data/JTM_upgrade/spatial_gse210616"
OUT  = ROOT / "02_results/JTM_upgrade/spatial_gse210616"
OUT.mkdir(parents=True, exist_ok=True)

fz = json.loads((ROOT / "02_results/JTM_upgrade/model_freeze/frozen_coefficients.json").read_text())
W_RIDGE = fz["weights"]
GMSIG4_GENES = ["IFNG", "IDO1", "CXCL10", "ENO1"]

# Strictly disjoint GMSig-4-free gene sets (Protocol Stage C)
EPI_GENES = ["EPCAM", "KRT8", "KRT18", "KRT19", "CDH1"]
STR_GENES = ["COL1A1", "COL1A2", "COL3A1", "DCN", "LUM"]
CD8_GENES = ["CD3D", "CD3E", "CD8A", "CD8B"]
ALL_NEEDED = set(GMSIG4_GENES + EPI_GENES + STR_GENES + CD8_GENES + ["STAT1", "GZMA", "PRF1", "CXCL12", "CXCR4"])


def load_visium_h5(h5_path: Path):
    """Load 10x filtered_feature_bc_matrix.h5 and return (barcodes, gene_names, csc_matrix [genes x spots])."""
    with h5py.File(h5_path, "r") as f:
        grp = f["matrix"]
        data = grp["data"][:]
        indices = grp["indices"][:]
        indptr = grp["indptr"][:]
        shape = tuple(grp["shape"][:])
        barcodes = [b.decode("utf-8") if isinstance(b, bytes) else str(b) for b in grp["barcodes"][:]]
        feat = grp["features"]
        names = [g.decode("utf-8") if isinstance(g, bytes) else str(g) for g in feat["name"][:]]
    mat = csc_matrix((data, indices, indptr), shape=shape)
    return barcodes, names, mat


def load_positions(pos_gz: Path) -> pd.DataFrame:
    with gzip.open(pos_gz, "rt", encoding="utf-8", errors="ignore") as f:
        first = f.readline().strip()
        f.seek(0)
        if first.lower().startswith("barcode"):
            df = pd.read_csv(f)
        else:
            df = pd.read_csv(f, header=None, names=["barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres"])
    return df


def visium_hex_margin(rows: np.ndarray, cols: np.ndarray, nest_mask: np.ndarray, max_steps: int = 2) -> np.ndarray:
    """Identify non-nest spots within <= max_steps Visium hexagonal lattice steps of a Tumor_Nest spot.
    In 10x Visium odd-r offset coordinates (array_row, array_col), each 1-step hex neighbor satisfies
    (|dr| == 0 & |dc| == 2) or (|dr| == 1 & |dc| == 1); 2-step hex distance is exact via BFS on the hex graph."""
    coord_to_idx = {(int(r), int(c)): i for i, (r, c) in enumerate(zip(rows, cols))}
    nbr_deltas = [(0, -2), (0, 2), (-1, -1), (-1, 1), (1, -1), (1, 1)]
    frontier = set(np.where(nest_mask)[0].tolist())
    visited = set(frontier)
    margin = set()
    for _ in range(max_steps):
        nxt = set()
        for idx in frontier:
            r, c = int(rows[idx]), int(cols[idx])
            for dr, dc in nbr_deltas:
                j = coord_to_idx.get((r + dr, c + dc))
                if j is not None and j not in visited:
                    visited.add(j)
                    nxt.add(j)
                    if not nest_mask[j]:
                        margin.add(j)
        frontier = nxt
    mask = np.zeros(len(rows), dtype=bool)
    if margin:
        mask[list(margin)] = True
    return mask


def parse_nac_table():
    """Extract patient-level NAC status (chemo-naive vs NAC) from can-22-2682 supplement."""
    xlsx = list(DATA.glob("can-22-2682*.xlsx"))[0]
    df = pd.read_excel(xlsx)
    # Forward-fill patient info only across consecutive non-blank rows? Better: keep rows with numeric Patient ID
    out = {}
    for _, row in df.iterrows():
        pid = row.get("Patient ID")
        try:
            pnum = int(float(pid))
        except Exception:
            continue
        tx = row.get("Neoadjuvant tx")
        is_nac = pd.notna(tx) and str(tx).strip() != "" and str(tx).strip().lower() not in ("none", "na", "nan")
        out[pnum] = {
            "nac_status": "NAC" if is_nac else "Chemo_Naive",
            "neoadjuvant_tx": "" if pd.isna(tx) else str(tx).strip(),
            "race": "" if pd.isna(row.get("Race")) else str(row.get("Race")).strip(),
            "stage": "" if pd.isna(row.get("Stage")) else str(row.get("Stage")).strip(),
            "age": pd.to_numeric(row.get("Age"), errors="coerce"),
        }
    return out


def main():
    meta = pd.read_csv(DATA / "gse210616_gsm_metadata.csv")
    nac_map = parse_nac_table()
    print(f"Loaded GSE210616 metadata: {len(meta)} sections, {meta['patient_id'].nunique()} unique patients; NAC parsed for {len(nac_map)} patients.")

    sec_rows = []
    spot_cache = {}  # gsm -> DataFrame for representative slide export

    for _, mrow in meta.iterrows():
        gsm = mrow["gsm"]
        sdir = DATA / "raw" / gsm
        h5_path = sdir / mrow["h5_file"]
        pos_path = sdir / mrow["positions_file"]
        if not h5_path.exists() or not pos_path.exists():
            print(f"WARNING: Missing files for {gsm}, skipping.")
            continue

        barcodes, gnames, mat = load_visium_h5(h5_path)
        pos_df = load_positions(pos_path)
        pos_df = pos_df[pos_df["in_tissue"] == 1].copy()
        bc_to_col = {bc: i for i, bc in enumerate(barcodes)}
        pos_df["col_idx"] = pos_df["barcode"].map(bc_to_col)
        pos_df = pos_df[pos_df["col_idx"].notna()].copy()
        pos_df["col_idx"] = pos_df["col_idx"].astype(int)
        col_idx = pos_df["col_idx"].values

        sub_mat = mat[:, col_idx]
        total_umi = np.asarray(sub_mat.sum(axis=0)).ravel().astype(float)
        valid_spots = total_umi > 0
        sub_mat = sub_mat[:, valid_spots]
        total_umi = total_umi[valid_spots]
        pos_df = pos_df.iloc[valid_spots].reset_index(drop=True)
        n_spots = len(pos_df)

        # Map needed genes (sum duplicates if any)
        g_idx_map = {}
        for gi, gn in enumerate(gnames):
            if gn in ALL_NEEDED:
                g_idx_map.setdefault(gn, []).append(gi)

        # Per-spot log1p(CP10K) for needed genes
        spot_expr = {}
        pb_counts = {}
        total_pb_umi = float(total_umi.sum())
        for g in ALL_NEEDED:
            if g in g_idx_map:
                raw_g = np.asarray(sub_mat[g_idx_map[g], :].sum(axis=0)).ravel().astype(float)
                spot_expr[g] = np.log1p(raw_g / total_umi * 1e4)
                pb_counts[g] = np.log1p(float(raw_g.sum()) / total_pb_umi * 1e4)
            else:
                spot_expr[g] = np.zeros(n_spots, dtype=float)
                pb_counts[g] = 0.0

        epi_s = np.mean([spot_expr[g] for g in EPI_GENES], axis=0)
        str_s = np.mean([spot_expr[g] for g in STR_GENES], axis=0)
        cd8_s = np.mean([spot_expr[g] for g in CD8_GENES], axis=0)
        d_epi_str = epi_s - str_s

        nest_mask = d_epi_str > 0
        margin_mask = visium_hex_margin(pos_df["array_row"].values, pos_df["array_col"].values, nest_mask, max_steps=2)

        c_core = float(np.mean(cd8_s[nest_mask])) if nest_mask.sum() > 0 else np.nan
        c_margin = float(np.mean(cd8_s[margin_mask])) if margin_mask.sum() > 0 else np.nan
        c_overall = float(np.mean(cd8_s))
        cd8_pi = float(np.log2((c_core + 0.01) / (c_margin + 0.01))) if (np.isfinite(c_core) and np.isfinite(c_margin)) else np.nan

        pnum = int(mrow["patient_num"])
        ninfo = nac_map.get(pnum, {"nac_status": "Unknown", "neoadjuvant_tx": "", "race": "", "stage": "", "age": np.nan})

        sec_rows.append({
            "gsm": gsm,
            "sample_code": mrow["sample_code"],
            "patient_num": pnum,
            "patient_id": mrow["patient_id"],
            "section_num": int(mrow["section_num"]),
            "cohort": mrow["cohort"],
            "nac_status": ninfo["nac_status"],
            "stage": ninfo["stage"],
            "n_spots": n_spots,
            "n_nest_spots": int(nest_mask.sum()),
            "n_margin_spots": int(margin_mask.sum()),
            "frac_nest": float(nest_mask.mean()),
            "C_core": c_core,
            "C_margin": c_margin,
            "C_overall": c_overall,
            "CD8_PI": cd8_pi,
            "pb_IFNG": pb_counts["IFNG"],
            "pb_IDO1": pb_counts["IDO1"],
            "pb_CXCL10": pb_counts["CXCL10"],
            "pb_ENO1": pb_counts["ENO1"],
            "pb_STAT1": pb_counts["STAT1"],
            "pb_GZMA": pb_counts["GZMA"],
            "pb_PRF1": pb_counts["PRF1"],
            "pb_CXCL12": pb_counts["CXCL12"],
            "mean_epi": float(np.mean(epi_s)),
            "mean_str": float(np.mean(str_s)),
        })

        # Keep lightweight spot data for representative slides
        spot_cache[gsm] = pd.DataFrame({
            "barcode": pos_df["barcode"].values,
            "array_row": pos_df["array_row"].values,
            "array_col": pos_df["array_col"].values,
            "pxl_row": pos_df["pxl_row_in_fullres"].values,
            "pxl_col": pos_df["pxl_col_in_fullres"].values,
            "compartment": np.where(nest_mask, "Tumor_Nest", np.where(margin_mask, "Peritumoral_Margin", "Distal_Stroma")),
            "epi_score": epi_s,
            "str_score": str_s,
            "cd8_score": cd8_s,
            "ENO1": spot_expr["ENO1"],
            "CXCL10": spot_expr["CXCL10"],
            "IDO1": spot_expr["IDO1"],
            "IFNG": spot_expr["IFNG"],
        })

    sec_df = pd.DataFrame(sec_rows)
    sec_df.to_csv(OUT / "section_level_spatial_metrics.csv", index=False)

    # Aggregate to patient level (N = 22 independent patients)
    num_cols = [
        "n_spots", "n_nest_spots", "n_margin_spots", "frac_nest",
        "C_core", "C_margin", "C_overall", "CD8_PI",
        "pb_IFNG", "pb_IDO1", "pb_CXCL10", "pb_ENO1", "pb_STAT1", "pb_GZMA", "pb_PRF1", "pb_CXCL12",
        "mean_epi", "mean_str",
    ]
    pt_df = sec_df.groupby(["patient_id", "patient_num", "cohort", "nac_status", "stage"], as_index=False).agg(
        n_sections=("gsm", "count"),
        **{c: (c, "mean") for c in num_cols}
    )
    # Compute frozen GMSig-4 across N=22 patients (gene-wise sample SD ddof=1)
    for g in GMSIG4_GENES:
        vals = pt_df[f"pb_{g}"].values
        pt_df[f"Z_{g}"] = (vals - np.mean(vals)) / np.std(vals, ddof=1)
    pt_df["GMSig4"] = sum(W_RIDGE[g] * pt_df[f"Z_{g}"] for g in GMSIG4_GENES)
    pt_df["Immune3_frozen"] = sum(W_RIDGE[g] * pt_df[f"Z_{g}"] for g in ["IFNG", "IDO1", "CXCL10"])
    pt_df["ENO1_only"] = pt_df["Z_ENO1"]
    pt_df["RiskGroup"] = np.where(pt_df["GMSig4"] > pt_df["GMSig4"].median(), "High", "Low")
    pt_df = pt_df.sort_values("GMSig4").reset_index(drop=True)
    pt_df.to_csv(OUT / "patient_level_spatial_metrics.csv", index=False)

    # Statistical summary across All (N=22) and Chemo-Naive sensitivity (N=13)
    stat_rows = []
    for sub_name, sub in [("All_TNBC_n22", pt_df), ("Chemo_Naive_n13", pt_df[pt_df["nac_status"] == "Chemo_Naive"].copy())]:
        if len(sub) < 5:
            continue
        sub_grp = np.where(sub["GMSig4"] > sub["GMSig4"].median(), "High", "Low")
        for score_col in ["GMSig4", "Immune3_frozen", "ENO1_only"]:
            for metric in ["CD8_PI", "C_overall", "C_core", "C_margin", "pb_STAT1", "pb_GZMA", "pb_CXCL12", "frac_nest"]:
                valid = sub[[score_col, metric]].dropna()
                rho, p_sp = spearmanr(valid[score_col], valid[metric])
                hi = sub.loc[sub_grp == "High", metric].dropna()
                lo = sub.loc[sub_grp == "Low", metric].dropna()
                _, p_w = mannwhitneyu(hi, lo, alternative="two-sided") if (len(hi) > 1 and len(lo) > 1) else (np.nan, np.nan)
                stat_rows.append({
                    "subset": sub_name, "n_patients": len(valid), "score": score_col, "metric": metric,
                    "spearman_rho": rho, "spearman_p": p_sp,
                    "median_Low": float(np.median(lo)), "median_High": float(np.median(hi)),
                    "wilcoxon_p": p_w,
                })
    stat_df = pd.DataFrame(stat_rows)
    stat_df.to_csv(OUT / "spatial_statistical_summary.csv", index=False)

    # Export spot tables for 1 representative Low-GMSig4 and 1 High-GMSig4 chemo-naive patient
    cn = pt_df[pt_df["nac_status"] == "Chemo_Naive"].sort_values("GMSig4")
    low_pt = cn.iloc[0]["patient_id"]
    high_pt = cn.iloc[-1]["patient_id"]
    rep_info = []
    for tag, pid in [("Low_GMSig4_Rep", low_pt), ("High_GMSig4_Rep", high_pt)]:
        sec_sub = sec_df[sec_df["patient_id"] == pid].sort_values("n_spots", ascending=False).iloc[0]
        gsm = sec_sub["gsm"]
        spot_cache[gsm].to_csv(OUT / f"spots_{tag}_{pid}_{gsm}.csv", index=False)
        rep_info.append({"role": tag, "patient_id": pid, "gsm": gsm, "sample_code": sec_sub["sample_code"],
                         "GMSig4": float(pt_df.loc[pt_df["patient_id"] == pid, "GMSig4"].values[0]),
                         "C_overall": float(sec_sub["C_overall"]), "C_core": float(sec_sub["C_core"]),
                         "C_margin": float(sec_sub["C_margin"]), "CD8_PI": float(sec_sub["CD8_PI"])})
    pd.DataFrame(rep_info).to_csv(OUT / "representative_slides.csv", index=False)

    print("\n=== GSE210616 Patient-Level Spatial Summary (N=22 & Chemo-Naive) ===")
    print(stat_df[(stat_df["score"] == "GMSig4") & (stat_df["metric"].isin(["CD8_PI", "C_overall", "C_core", "C_margin", "pb_STAT1", "pb_GZMA", "pb_CXCL12"]))].to_string(index=False))
    print("\nRepresentative slides:", rep_info)


if __name__ == "__main__":
    main()
