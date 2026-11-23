#!/usr/bin/env python3
"""
Stage C2: Narrow CPTAC Proteomics Orthogonal Validation (brca_cptac_2020; Krug et al., Cell 2020)
for the GMSig-4 JTM upgrade.

Strict methodological rules:
- Exact double-precision frozen weights from frozen_coefficients.json.
- data_mrna_seq_fpkm.txt is already log2-transformed (log2-FPKM, containing negative values);
  do NOT apply a second log2(x+1). Apply gene-wise cohort Z-scaling with sample SD (ddof=1).
- Median split: strictly greater than cohort median = High; otherwise Low.
- Zero imputation for mass-spec proteomics; NO 4-protein GMSig-4 score (IFNG = 0/28).
- Primary cohort: TNBC_UPDATED_CLINICAL_STATUS == "Positive" (n = 28).
- Sensitivity cohort: strict clinical IHC ER-/PR-/HER2- (n = 16).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from scipy import stats

BASE_DIR = Path(r"D:\Bioinformatics\TNBC_Microbiome")
DATA_DIR = BASE_DIR / "01_data" / "JTM_upgrade" / "cptac"
OUT_DIR = BASE_DIR / "02_results" / "JTM_upgrade" / "cptac"
MODEL_JSON = BASE_DIR / "02_results" / "JTM_upgrade" / "model_freeze" / "frozen_coefficients.json"

GMSIG4_GENES = ["IFNG", "IDO1", "CXCL10", "ENO1"]
GLYCOLYSIS_PROXY_GENES = ["LDHA", "PKM", "GAPDH", "PGK1"]
IMMUNE_PROXY_GENES = ["STAT1", "GZMA", "PRF1", "CD8A"]
ALL_CHECK_GENES = ["IFNG", "IDO1", "CXCL10", "ENO1", "STAT1", "CD274"] + GLYCOLYSIS_PROXY_GENES + ["GZMA", "PRF1", "CD8A"]


def zscore_series(s: pd.Series) -> pd.Series:
    vals = pd.to_numeric(s, errors="coerce")
    mean = vals.mean(skipna=True)
    sd = vals.std(ddof=1, skipna=True)
    return (vals - mean) / sd


def compute_gmsig4(df_rna: pd.DataFrame, weights: dict[str, float]) -> tuple[pd.Series, pd.DataFrame]:
    z_df = pd.DataFrame(index=df_rna.index)
    score = pd.Series(0.0, index=df_rna.index, dtype=float)
    for g in GMSIG4_GENES:
        z = zscore_series(df_rna[f"RNA_{g}"])
        z_df[f"RNA_Z_{g}"] = z
        score = score + weights[g] * z
    return score, z_df


def spearman_and_wilcoxon(
    df: pd.DataFrame,
    score_col: str,
    group_col: str,
    target_col: str,
    cohort_label: str,
    category: str,
    feature_label: str,
) -> dict[str, Any]:
    sub = df[[score_col, group_col, target_col]].dropna()
    n = len(sub)
    if n < 3:
        return {
            "cohort": cohort_label,
            "analysis_type": "GMSig4_association",
            "category": category,
            "feature": feature_label,
            "n": n,
            "n_high": int((sub[group_col] == "High").sum()) if n else 0,
            "n_low": int((sub[group_col] == "Low").sum()) if n else 0,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "median_high": np.nan,
            "median_low": np.nan,
            "median_diff_high_minus_low": np.nan,
            "wilcoxon_u": np.nan,
            "wilcoxon_p": np.nan,
            "note": "Insufficient non-missing observations",
        }

    rho, sp_p = stats.spearmanr(sub[score_col], sub[target_col])
    r, pr_p = stats.pearsonr(sub[score_col], sub[target_col])

    hi = sub.loc[sub[group_col] == "High", target_col]
    lo = sub.loc[sub[group_col] == "Low", target_col]
    if len(hi) >= 1 and len(lo) >= 1:
        u_stat, wx_p = stats.mannwhitneyu(hi, lo, alternative="two-sided")
        med_hi = float(hi.median())
        med_lo = float(lo.median())
        med_diff = med_hi - med_lo
    else:
        u_stat, wx_p, med_hi, med_lo, med_diff = np.nan, np.nan, np.nan, np.nan, np.nan

    return {
        "cohort": cohort_label,
        "analysis_type": "GMSig4_association",
        "category": category,
        "feature": feature_label,
        "n": int(n),
        "n_high": int(len(hi)),
        "n_low": int(len(lo)),
        "spearman_rho": float(rho),
        "spearman_p": float(sp_p),
        "pearson_r": float(r),
        "pearson_p": float(pr_p),
        "median_high": med_hi,
        "median_low": med_lo,
        "median_diff_high_minus_low": med_diff,
        "wilcoxon_u": float(u_stat) if not np.isnan(u_stat) else np.nan,
        "wilcoxon_p": float(wx_p) if not np.isnan(wx_p) else np.nan,
        "note": "",
    }


def rna_protein_concordance(
    df: pd.DataFrame, gene: str, cohort_label: str, total_n: int
) -> dict[str, Any]:
    r_col = f"RNA_{gene}"
    p_col = f"PROT_{gene}"
    if p_col not in df.columns or df[p_col].notna().sum() == 0:
        return {
            "cohort": cohort_label,
            "analysis_type": "RNA_protein_concordance",
            "category": "RNA_to_Protein",
            "feature": gene,
            "n": 0,
            "n_high": np.nan,
            "n_low": np.nan,
            "spearman_rho": np.nan,
            "spearman_p": np.nan,
            "pearson_r": np.nan,
            "pearson_p": np.nan,
            "median_high": np.nan,
            "median_low": np.nan,
            "median_diff_high_minus_low": np.nan,
            "wilcoxon_u": np.nan,
            "wilcoxon_p": np.nan,
            "note": f"Not detected in MS proteomics (0/{total_n})",
        }

    sub = df[[r_col, p_col]].dropna()
    n = len(sub)
    rho, sp_p = stats.spearmanr(sub[r_col], sub[p_col])
    r, pr_p = stats.pearsonr(sub[r_col], sub[p_col])
    n_prot = int(df[p_col].notna().sum())
    n_rna = int(df[r_col].notna().sum())
    return {
        "cohort": cohort_label,
        "analysis_type": "RNA_protein_concordance",
        "category": "RNA_to_Protein",
        "feature": gene,
        "n": int(n),
        "n_high": np.nan,
        "n_low": np.nan,
        "spearman_rho": float(rho),
        "spearman_p": float(sp_p),
        "pearson_r": float(r),
        "pearson_p": float(pr_p),
        "median_high": np.nan,
        "median_low": np.nan,
        "median_diff_high_minus_low": np.nan,
        "wilcoxon_u": np.nan,
        "wilcoxon_p": np.nan,
        "note": f"RNA detected {n_rna}/{total_n}; Protein detected {n_prot}/{total_n}; non-imputed complete pairs {n}/{total_n}",
    }


def mnar_missingness_check(df: pd.DataFrame, gene: str, cohort_label: str) -> dict[str, Any] | None:
    r_col = f"RNA_{gene}"
    p_col = f"PROT_{gene}"
    if p_col not in df.columns:
        return None
    det = df.loc[df[p_col].notna(), r_col].dropna()
    mis = df.loc[df[p_col].isna(), r_col].dropna()
    if len(det) == 0 or len(mis) == 0:
        return None
    u_stat, p_val = stats.mannwhitneyu(det, mis, alternative="two-sided")
    return {
        "cohort": cohort_label,
        "analysis_type": "MNAR_missingness_check",
        "category": "Protein_Missingness_vs_RNA",
        "feature": f"{gene}_RNA_in_ProtDetected_vs_ProtNA",
        "n": int(len(det) + len(mis)),
        "n_high": int(len(det)),
        "n_low": int(len(mis)),
        "spearman_rho": np.nan,
        "spearman_p": np.nan,
        "pearson_r": np.nan,
        "pearson_p": np.nan,
        "median_high": float(det.median()),
        "median_low": float(mis.median()),
        "median_diff_high_minus_low": float(det.median() - mis.median()),
        "wilcoxon_u": float(u_stat),
        "wilcoxon_p": float(p_val),
        "note": "n_high/median_high = Protein Detected; n_low/median_low = Protein NA (tests abundance-driven MNAR)",
    }


def main() -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)

    with MODEL_JSON.open("r", encoding="utf-8") as f:
        frozen = json.load(f)
    weights: dict[str, float] = frozen["weights"]
    print("Loaded frozen weights:", weights)

    # 1. Load clinical tables
    pat = pd.read_csv(DATA_DIR / "data_clinical_patient.txt", sep="\t", comment="#")
    smp = pd.read_csv(DATA_DIR / "data_clinical_sample.txt", sep="\t", comment="#")
    clin = pat.merge(smp, on="PATIENT_ID", how="inner")
    assert len(clin) == 122, f"Expected 122 CPTAC patients, got {len(clin)}"

    # Confirm 0 overlap with TCGA-BRCA IDs
    tcga_overlap = clin["PATIENT_ID"].astype(str).str.startswith("TCGA").sum()
    assert tcga_overlap == 0, "Unexpected TCGA overlap!"

    # Subset confirmed TNBC (n=28) and flag strict clinical IHC subset (n=16)
    tnbc = clin[clin["TNBC_UPDATED_CLINICAL_STATUS"] == "Positive"].copy().reset_index(drop=True)
    assert len(tnbc) == 28, f"Expected 28 TNBC patients, got {len(tnbc)}"

    tnbc["STRICT_IHC_TNBC"] = (
        (tnbc["ER_UPDATED_CLINICAL_STATUS"].astype(str).str.lower() == "negative")
        & (tnbc["PR_CLINICAL_STATUS"].astype(str).str.lower() == "negative")
        & (tnbc["ERBB2_UPDATED_CLINICAL_STATUS"].astype(str).str.lower() == "negative")
    )
    assert int(tnbc["STRICT_IHC_TNBC"].sum()) == 16, f"Expected 16 strict IHC TNBC, got {tnbc['STRICT_IHC_TNBC'].sum()}"

    # Parse numeric CD3_TILS_COUNTS
    tnbc["CD3_TILS_COUNTS_NUM"] = pd.to_numeric(tnbc["CD3_TILS_COUNTS"], errors="coerce")
    assert int(tnbc["CD3_TILS_COUNTS_NUM"].notna().sum()) == 13, "Expected 13 numeric CD3_TILS_COUNTS"

    sample_ids = tnbc["SAMPLE_ID"].tolist()

    # 2. Load RNA-seq (already log2-transformed FPKM)
    rna = pd.read_csv(DATA_DIR / "data_mrna_seq_fpkm.txt", sep="\t")
    rna_sym_col = rna.columns[0]
    rna_sub = rna[rna[rna_sym_col].isin(ALL_CHECK_GENES)].set_index(rna_sym_col)
    for g in ALL_CHECK_GENES:
        assert g in rna_sub.index, f"Gene {g} missing from RNA matrix"
        tnbc[f"RNA_{g}"] = [pd.to_numeric(rna_sub.loc[g, sid], errors="coerce") for sid in sample_ids]

    # 3. Load mass-spec proteomics (TMT-10-plex log2 ratios, non-imputed)
    prot = pd.read_csv(DATA_DIR / "data_protein_quantification.txt", sep="\t")
    prot_ref_col = prot.columns[0]
    prot["GENE"] = prot[prot_ref_col].astype(str).str.split("|").str[0]

    # Explicitly confirm IFNG = 0/122 and 0/28, CD274 = 0/122 and 0/28
    assert (prot["GENE"] == "IFNG").sum() == 0, "IFNG unexpectedly found in proteomics!"
    assert (prot["GENE"] == "CD274").sum() == 0, "CD274 unexpectedly found in proteomics!"

    prot_ref_map: dict[str, str] = {}
    for g in ALL_CHECK_GENES:
        rows_g = prot[prot["GENE"] == g]
        if len(rows_g) == 0:
            tnbc[f"PROT_{g}"] = np.nan
            prot_ref_map[g] = "ABSENT (0/122 BRCA, 0/28 TNBC)"
        else:
            assert len(rows_g) == 1, f"Multiple protein rows for {g}: {rows_g[prot_ref_col].tolist()}"
            row_g = rows_g.iloc[0]
            prot_ref_map[g] = str(row_g[prot_ref_col])
            tnbc[f"PROT_{g}"] = [pd.to_numeric(row_g[sid], errors="coerce") for sid in sample_ids]

    # Confirm exact protein detection counts in TNBC
    assert int(tnbc["PROT_IFNG"].notna().sum()) == 0
    assert int(tnbc["PROT_ENO1"].notna().sum()) == 28
    assert int(tnbc["PROT_IDO1"].notna().sum()) == 23
    assert int(tnbc["PROT_CXCL10"].notna().sum()) == 21
    assert int(tnbc["PROT_STAT1"].notna().sum()) == 28
    assert int(tnbc["PROT_CD274"].notna().sum()) == 0

    # 4. Compute frozen RNA GMSig-4 score across n=28 TNBC
    tnbc["GMSig4_Score"], z_df = compute_gmsig4(tnbc, weights)
    for col in z_df.columns:
        tnbc[col] = z_df[col]
    med_28 = float(tnbc["GMSig4_Score"].median())
    tnbc["GMSig4_Group"] = np.where(tnbc["GMSig4_Score"] > med_28, "High", "Low")
    assert int((tnbc["GMSig4_Group"] == "High").sum()) == 14
    assert int((tnbc["GMSig4_Group"] == "Low").sum()) == 14

    # Opposing arms for biological transparency (non-circular decomposition)
    tnbc["GMSig4_Immune3_Arm"] = (
        weights["IFNG"] * tnbc["RNA_Z_IFNG"]
        + weights["IDO1"] * tnbc["RNA_Z_IDO1"]
        + weights["CXCL10"] * tnbc["RNA_Z_CXCL10"]
    )
    tnbc["GMSig4_ENO1_Arm"] = weights["ENO1"] * tnbc["RNA_Z_ENO1"]

    # 5. Compute non-imputed orthogonal protein proxies (within n=28 TNBC)
    # Glycolysis hallmark protein proxy (strictly excluding ENO1): mean Z of LDHA, PKM, GAPDH, PGK1 (all 28/28)
    glyc_z_cols = []
    for g in GLYCOLYSIS_PROXY_GENES:
        col_z = f"PROT_Z_{g}"
        tnbc[col_z] = zscore_series(tnbc[f"PROT_{g}"])
        glyc_z_cols.append(col_z)
    tnbc["PROT_Glycolysis_Proxy_MeanZ"] = tnbc[glyc_z_cols].mean(axis=1)

    # Immune protein proxy:
    # (1) Complete-case 2-protein (STAT1 + GZMA, both 28/28)
    # (2) Available-case 4-protein (mean Z across non-missing of STAT1, GZMA, PRF1, CD8A; n=28)
    # (3) Strict complete-case 4-protein (STAT1, GZMA, PRF1, CD8A all non-missing)
    imm_z_cols = []
    for g in IMMUNE_PROXY_GENES:
        col_z = f"PROT_Z_{g}"
        tnbc[col_z] = zscore_series(tnbc[f"PROT_{g}"])
        imm_z_cols.append(col_z)
    tnbc["PROT_Immune_Proxy_STAT1_GZMA_MeanZ"] = tnbc[["PROT_Z_STAT1", "PROT_Z_GZMA"]].mean(axis=1)
    tnbc["PROT_Immune_Proxy_4Avail_MeanZ"] = tnbc[imm_z_cols].mean(axis=1, skipna=True)
    tnbc["PROT_Immune_Proxy_4Complete_MeanZ"] = tnbc[imm_z_cols].mean(axis=1, skipna=False)

    # 6. Build summary results table across:
    #    (i) TNBC_All (n=28)
    #    (ii) TNBC_Strict_IHC_Frozen28Score (n=16, using the 28-TNBC frozen score & median split)
    #    (iii) TNBC_Strict_IHC_ReZ16 (n=16, re-Z-scaled within n=16 sensitivity)
    summary_rows: list[dict[str, Any]] = []

    strict_df = tnbc[tnbc["STRICT_IHC_TNBC"]].copy().reset_index(drop=True)
    strict_df["GMSig4_Score_ReZ16"], _ = compute_gmsig4(strict_df, weights)
    med_16 = float(strict_df["GMSig4_Score_ReZ16"].median())
    strict_df["GMSig4_Group_ReZ16"] = np.where(strict_df["GMSig4_Score_ReZ16"] > med_16, "High", "Low")
    # Recompute protein proxies Z within strict n=16 for ReZ16 sensitivity
    strict_df["PROT_Glycolysis_Proxy_MeanZ_16"] = pd.concat(
        [zscore_series(strict_df[f"PROT_{g}"]) for g in GLYCOLYSIS_PROXY_GENES], axis=1
    ).mean(axis=1)
    strict_df["PROT_Immune_Proxy_STAT1_GZMA_MeanZ_16"] = pd.concat(
        [zscore_series(strict_df[f"PROT_{g}"]) for g in ["STAT1", "GZMA"]], axis=1
    ).mean(axis=1)
    strict_df["PROT_Immune_Proxy_4Avail_MeanZ_16"] = pd.concat(
        [zscore_series(strict_df[f"PROT_{g}"]) for g in IMMUNE_PROXY_GENES], axis=1
    ).mean(axis=1, skipna=True)

    # A. RNA-to-protein concordance (n=28 and n=16)
    for cohort_label, df_c in [("TNBC_All_n28", tnbc), ("TNBC_StrictIHC_n16", strict_df)]:
        for g in ["ENO1", "IDO1", "CXCL10", "IFNG", "STAT1", "CD274", "LDHA", "PKM", "GAPDH", "PGK1", "GZMA", "PRF1", "CD8A"]:
            summary_rows.append(rna_protein_concordance(df_c, g, cohort_label, len(df_c)))

    # B. MNAR missingness check for IDO1 and CXCL10 (and PRF1)
    for g in ["IDO1", "CXCL10", "PRF1"]:
        res = mnar_missingness_check(tnbc, g, "TNBC_All_n28")
        if res is not None:
            summary_rows.append(res)

    # C. GMSig-4 associations (n=28 primary)
    targets_28 = [
        ("GMSig4_Proteins", "PROT_ENO1", "ENO1_protein"),
        ("GMSig4_Proteins", "PROT_IDO1", "IDO1_protein"),
        ("GMSig4_Proteins", "PROT_CXCL10", "CXCL10_protein"),
        ("Clinical_TME", "ESTIMATE_IMMUNE_SCORE", "ESTIMATE_IMMUNE_SCORE"),
        ("Clinical_TME", "ESTIMATE_STROMAL_SCORE", "ESTIMATE_STROMAL_SCORE"),
        ("Clinical_TME", "ESTIMATE_TUMORPURITY", "ESTIMATE_TUMORPURITY"),
        ("Clinical_TME", "CIBERSORT_ABSOLUTE_SCORE", "CIBERSORT_ABSOLUTE_SCORE"),
        ("Clinical_TME", "XCELL_IMMUNE_SCORE", "XCELL_IMMUNE_SCORE"),
        ("Pathology_TILs", "CD3_TILS_COUNTS_NUM", "CD3_TILS_COUNTS"),
        ("Glycolysis_Proteins", "PROT_LDHA", "LDHA_protein"),
        ("Glycolysis_Proteins", "PROT_PKM", "PKM_protein"),
        ("Glycolysis_Proteins", "PROT_GAPDH", "GAPDH_protein"),
        ("Glycolysis_Proteins", "PROT_PGK1", "PGK1_protein"),
        ("Glycolysis_Proteins", "PROT_Glycolysis_Proxy_MeanZ", "Glycolysis_4Prot_Proxy_MeanZ_noENO1"),
        ("Immune_Proteins", "PROT_STAT1", "STAT1_protein"),
        ("Immune_Proteins", "PROT_GZMA", "GZMA_protein"),
        ("Immune_Proteins", "PROT_PRF1", "PRF1_protein"),
        ("Immune_Proteins", "PROT_CD8A", "CD8A_protein"),
        ("Immune_Proteins", "PROT_Immune_Proxy_STAT1_GZMA_MeanZ", "Immune_2Prot_Complete_STAT1_GZMA_MeanZ"),
        ("Immune_Proteins", "PROT_Immune_Proxy_4Avail_MeanZ", "Immune_4Prot_Avail_MeanZ"),
        ("Immune_Proteins", "PROT_Immune_Proxy_4Complete_MeanZ", "Immune_4Prot_Complete_MeanZ"),
    ]
    for cat, col, label in targets_28:
        summary_rows.append(
            spearman_and_wilcoxon(tnbc, "GMSig4_Score", "GMSig4_Group", col, "TNBC_All_n28", cat, label)
        )

    # D. Strict IHC TNBC sensitivity (n=16, both re-Z-scaled within n=16 and using full 28-TNBC score)
    targets_16 = [
        ("GMSig4_Proteins", "PROT_ENO1", "ENO1_protein"),
        ("GMSig4_Proteins", "PROT_IDO1", "IDO1_protein"),
        ("GMSig4_Proteins", "PROT_CXCL10", "CXCL10_protein"),
        ("Clinical_TME", "ESTIMATE_IMMUNE_SCORE", "ESTIMATE_IMMUNE_SCORE"),
        ("Pathology_TILs", "CD3_TILS_COUNTS_NUM", "CD3_TILS_COUNTS"),
        ("Glycolysis_Proteins", "PROT_Glycolysis_Proxy_MeanZ_16", "Glycolysis_4Prot_Proxy_MeanZ_noENO1"),
        ("Immune_Proteins", "PROT_STAT1", "STAT1_protein"),
        ("Immune_Proteins", "PROT_GZMA", "GZMA_protein"),
        ("Immune_Proteins", "PROT_Immune_Proxy_STAT1_GZMA_MeanZ_16", "Immune_2Prot_Complete_STAT1_GZMA_MeanZ"),
        ("Immune_Proteins", "PROT_Immune_Proxy_4Avail_MeanZ_16", "Immune_4Prot_Avail_MeanZ"),
    ]
    for cat, col, label in targets_16:
        summary_rows.append(
            spearman_and_wilcoxon(
                strict_df,
                "GMSig4_Score_ReZ16",
                "GMSig4_Group_ReZ16",
                col,
                "TNBC_StrictIHC_n16_ReZ16",
                cat,
                label,
            )
        )

    # Save patient table
    keep_cols = [
        "PATIENT_ID",
        "SAMPLE_ID",
        "AGE",
        "TUMOR_STAGE",
        "PAM50",
        "NMF_CLUSTER",
        "ER_UPDATED_CLINICAL_STATUS",
        "PR_CLINICAL_STATUS",
        "ERBB2_UPDATED_CLINICAL_STATUS",
        "ERBB2_PROTEOGENOMIC_STATUS",
        "TNBC_UPDATED_CLINICAL_STATUS",
        "STRICT_IHC_TNBC",
        "ESTIMATE_IMMUNE_SCORE",
        "ESTIMATE_STROMAL_SCORE",
        "ESTIMATE_TUMORPURITY",
        "CIBERSORT_ABSOLUTE_SCORE",
        "XCELL_IMMUNE_SCORE",
        "CD3_TILS_STATUS",
        "CD3_TILS_COUNTS",
        "CD3_TILS_COUNTS_NUM",
        "GMSig4_Score",
        "GMSig4_Group",
        "GMSig4_Immune3_Arm",
        "GMSig4_ENO1_Arm",
        "RNA_IFNG",
        "RNA_IDO1",
        "RNA_CXCL10",
        "RNA_ENO1",
        "RNA_Z_IFNG",
        "RNA_Z_IDO1",
        "RNA_Z_CXCL10",
        "RNA_Z_ENO1",
        "RNA_STAT1",
        "RNA_CD274",
        "PROT_IFNG",
        "PROT_IDO1",
        "PROT_CXCL10",
        "PROT_ENO1",
        "PROT_STAT1",
        "PROT_CD274",
        "PROT_LDHA",
        "PROT_PKM",
        "PROT_GAPDH",
        "PROT_PGK1",
        "PROT_GZMA",
        "PROT_PRF1",
        "PROT_CD8A",
        "PROT_Glycolysis_Proxy_MeanZ",
        "PROT_Immune_Proxy_STAT1_GZMA_MeanZ",
        "PROT_Immune_Proxy_4Avail_MeanZ",
        "PROT_Immune_Proxy_4Complete_MeanZ",
    ]
    patient_out = OUT_DIR / "cptac_tnbc_patient_table.csv"
    tnbc[keep_cols].to_csv(patient_out, index=False)

    summary_df = pd.DataFrame(summary_rows)
    summary_out = OUT_DIR / "cptac_protein_concordance_summary.csv"
    summary_df.to_csv(summary_out, index=False)

    # Print formatted results to stdout for verification
    print("\n=== COHORT SUMMARY ===")
    print(f"Total CPTAC BRCA: {len(clin)} | Confirmed TNBC: {len(tnbc)} | Strict IHC TNBC: {int(tnbc['STRICT_IHC_TNBC'].sum())}")
    print(f"PAM50 in TNBC: {tnbc['PAM50'].value_counts().to_dict()}")
    print(f"GMSig4_Score (n=28): min={tnbc['GMSig4_Score'].min():.4f}, median={med_28:.4f}, max={tnbc['GMSig4_Score'].max():.4f}, SD={tnbc['GMSig4_Score'].std(ddof=1):.4f}")
    print(f"GMSig4_Group (n=28): {tnbc['GMSig4_Group'].value_counts().to_dict()}")
    print(f"CD3_TILS_COUNTS numeric: n={int(tnbc['CD3_TILS_COUNTS_NUM'].notna().sum())} (High={int((tnbc.loc[tnbc['CD3_TILS_COUNTS_NUM'].notna(), 'GMSig4_Group']=='High').sum())}, Low={int((tnbc.loc[tnbc['CD3_TILS_COUNTS_NUM'].notna(), 'GMSig4_Group']=='Low').sum())})")
    print(f"Strict IHC CD3_TILS_COUNTS numeric: n={int(strict_df['CD3_TILS_COUNTS_NUM'].notna().sum())}")

    print("\n=== RNA-TO-PROTEIN CONCORDANCE (TNBC n=28) ===")
    r2p = summary_df[(summary_df["cohort"] == "TNBC_All_n28") & (summary_df["analysis_type"] == "RNA_protein_concordance")]
    print(r2p[["feature", "n", "spearman_rho", "spearman_p", "pearson_r", "pearson_p", "note"]].to_string(index=False))

    print("\n=== MNAR MISSINGNESS CHECK (TNBC n=28) ===")
    mnar = summary_df[summary_df["analysis_type"] == "MNAR_missingness_check"]
    print(mnar[["feature", "n_high", "n_low", "median_high", "median_low", "wilcoxon_p"]].to_string(index=False))

    print("\n=== GMSig-4 ASSOCIATIONS (TNBC n=28) ===")
    ga28 = summary_df[(summary_df["cohort"] == "TNBC_All_n28") & (summary_df["analysis_type"] == "GMSig4_association")]
    print(
        ga28[
            ["category", "feature", "n", "n_high", "n_low", "spearman_rho", "spearman_p", "median_high", "median_low", "wilcoxon_p"]
        ].to_string(index=False)
    )

    print("\n=== GMSig-4 ASSOCIATIONS (STRICT IHC TNBC n=16, ReZ16) ===")
    ga16 = summary_df[(summary_df["cohort"] == "TNBC_StrictIHC_n16_ReZ16") & (summary_df["analysis_type"] == "GMSig4_association")]
    print(
        ga16[
            ["category", "feature", "n", "n_high", "n_low", "spearman_rho", "spearman_p", "median_high", "median_low", "wilcoxon_p"]
        ].to_string(index=False)
    )

    print("\n=== RNA-TO-PROTEIN CONCORDANCE (STRICT IHC TNBC n=16) ===")
    r2p16 = summary_df[(summary_df["cohort"] == "TNBC_StrictIHC_n16") & (summary_df["analysis_type"] == "RNA_protein_concordance")]
    print(r2p16[r2p16["feature"].isin(["ENO1", "IDO1", "CXCL10", "STAT1", "IFNG", "CD274"])][["feature", "n", "spearman_rho", "spearman_p"]].to_string(index=False))

    print(f"\nSaved: {patient_out}")
    print(f"Saved: {summary_out}")


if __name__ == "__main__":
    main()
