import json
import os
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

hgmt_dir = os.path.join(project_root, "01_data", "HGMT")
res_dir  = os.path.join(project_root, "02_results")
os.makedirs(res_dir, exist_ok=True)

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranks = np.arange(1, n + 1)
    q = np.empty(n)
    q[order] = p[order] * n / ranks
    # Enforce monotonicity
    for i in range(n - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.clip(q, 0, 1)

# 1. Process 5 Case-Control cohorts in HGMT (Breast Neoplasms vs Health)
cc_projs = ["PRJNA658160", "PRJNA383849", "PRJNA453965", "PRJNA345373", "PRJNA751260"]
all_genus_diff = []
cohort_summary = []

for pid in cc_projs:
    runs = json.load(open(os.path.join(hgmt_dir, f"{pid}_runs.json"), encoding="utf-8"))
    det  = json.load(open(os.path.join(hgmt_dir, f"{pid}_details.json"), encoding="utf-8"))
    bact = json.load(open(os.path.join(hgmt_dir, f"{pid}_feature_bacteria.json"), encoding="utf-8"))

    run_pheno = {r["run_id"]: r.get("phenotype_name", "") for r in runs}
    # Keep Breast Neoplasms & Health
    keep_runs = {rid: ph for rid, ph in run_pheno.items() if ph in ("Breast Neoplasms", "Health")}

    # Extract genus-level relative abundance per run
    genus_map = defaultdict(lambda: defaultdict(float))
    for item in bact:
        rid = item["run_id"]
        if rid not in keep_runs:
            continue
        taxa = item["taxa"]
        # Parse last g__ token if present
        parts = taxa.split("|")
        g_part = [x for x in parts if x.startswith("g__")]
        if not g_part:
            continue
        gname = g_part[-1].replace("g__", "").strip()
        if not gname or gname in ("unclassified", "NA"):
            continue
        genus_map[rid][gname] += float(item["abundance"])

    valid_rids = [rid for rid in keep_runs if rid in genus_map]
    bc_rids = [r for r in valid_rids if keep_runs[r] == "Breast Neoplasms"]
    hc_rids = [r for r in valid_rids if keep_runs[r] == "Health"]
    all_genera = sorted({g for rid in valid_rids for g in genus_map[rid]})

    mat = pd.DataFrame(0.0, index=valid_rids, columns=all_genera)
    for rid in valid_rids:
        for g, val in genus_map[rid].items():
            mat.at[rid, g] = val

    # Save clean genus abundance matrix for this HGMT cohort
    mat_out = mat.copy()
    mat_out.insert(0, "phenotype", [keep_runs[r] for r in valid_rids])
    mat_out.insert(0, "project_id", pid)
    mat_out.index.name = "run_id"
    mat_out.to_csv(os.path.join(hgmt_dir, f"HGMT_genus_{pid}.csv"))

    # Alpha diversity (Shannon & Richness)
    vals = mat.values
    row_sums = vals.sum(axis=1, keepdims=True)
    probs = np.divide(vals, np.where(row_sums > 0, row_sums, 1.0))
    with np.errstate(divide="ignore", invalid="ignore"):
        shannon = -np.nansum(np.where(probs > 0, probs * np.log(probs), 0.0), axis=1)
    richness = (vals > 0).sum(axis=1)
    idx_bc = np.array([keep_runs[r] == "Breast Neoplasms" for r in valid_rids])
    idx_hc = ~idx_bc

    _, p_shan = stats.mannwhitneyu(shannon[idx_bc], shannon[idx_hc], alternative="two-sided")
    _, p_rich = stats.mannwhitneyu(richness[idx_bc], richness[idx_hc], alternative="two-sided")

    cohort_summary.append({
        "project_id": pid,
        "assay_type": det.get("assay_type"),
        "pmid": det.get("public_accession"),
        "n_BC": int(idx_bc.sum()),
        "n_Health": int(idx_hc.sum()),
        "n_genera": len(all_genera),
        "shannon_BC_med": float(np.median(shannon[idx_bc])),
        "shannon_Health_med": float(np.median(shannon[idx_hc])),
        "shannon_p": float(p_shan),
        "richness_BC_med": float(np.median(richness[idx_bc])),
        "richness_Health_med": float(np.median(richness[idx_hc])),
        "richness_p": float(p_rich),
    })

    # Differential genus analysis (prevalence >= 10%)
    rows_diff = []
    for g in all_genera:
        v = mat[g].values
        if (v > 0).mean() < 0.10:
            continue
        v_bc = v[idx_bc]
        v_hc = v[idx_hc]
        _, pval = stats.mannwhitneyu(v_bc, v_hc, alternative="two-sided")
        rows_diff.append({
            "project_id": pid,
            "assay_type": det.get("assay_type"),
            "genus": g,
            "mean_BC": float(v_bc.mean()),
            "mean_Health": float(v_hc.mean()),
            "prev_BC": float((v_bc > 0).mean()),
            "prev_Health": float((v_hc > 0).mean()),
            "log2FC": float(np.log2((v_bc.mean() + 1e-4) / (v_hc.mean() + 1e-4))),
            "wilcox_p": float(pval),
        })
    df_d = pd.DataFrame(rows_diff)
    if len(df_d) > 0:
        df_d["wilcox_fdr"] = bh_fdr(df_d["wilcox_p"].values)
        all_genus_diff.append(df_d)

df_csum = pd.DataFrame(cohort_summary)
df_csum.to_csv(os.path.join(res_dir, "HGMT_5cohorts_alpha_summary.csv"), index=False)
print("=== HGMT 5 Case-Control Breast Cancer Cohorts Summary ===")
print(df_csum.to_string(index=False))

df_gdiff = pd.concat(all_genus_diff, ignore_index=True)
df_gdiff.to_csv(os.path.join(res_dir, "HGMT_5cohorts_genus_diff.csv"), index=False)

# Cross-cohort consistency in HGMT
cons = df_gdiff.groupby("genus").agg(
    n_cohorts=("project_id", "count"),
    n_p05=("wilcox_p", lambda x: int((x < 0.05).sum())),
    n_fdr05=("wilcox_fdr", lambda x: int((x < 0.05).sum())),
    mean_log2FC=("log2FC", "mean"),
    wgs_p=("wilcox_p", lambda x: float(x[df_gdiff.loc[x.index, "assay_type"] == "WGS"].min()) if (df_gdiff.loc[x.index, "assay_type"] == "WGS").any() else np.nan),
    wgs_log2FC=("log2FC", lambda x: float(x[df_gdiff.loc[x.index, "assay_type"] == "WGS"].mean()) if (df_gdiff.loc[x.index, "assay_type"] == "WGS").any() else np.nan),
    cohorts_sig=("project_id", lambda x: ";".join(x[df_gdiff.loc[x.index, "wilcox_p"] < 0.05])),
).reset_index()
cons = cons.sort_values(["n_fdr05", "n_p05", "wgs_p"], ascending=[False, False, True])
cons.to_csv(os.path.join(res_dir, "HGMT_genus_cross_cohort_consistency.csv"), index=False)
print("\n=== Top 25 HGMT Cross-Cohort Differential Genera ===")
print(cons.head(25).to_string(index=False))

# 2. WGS KEGG KO differential analysis in PRJNA453965 (62 BC vs 71 Health)
print("\n=== 2. HGMT PRJNA453965 WGS KEGG KO Functional Analysis ===")
runs_wgs = json.load(open(os.path.join(hgmt_dir, "PRJNA453965_runs.json"), encoding="utf-8"))
ko_json  = json.load(open(os.path.join(hgmt_dir, "PRJNA453965_feature_ko.json"), encoding="utf-8"))
run_ph   = {r["run_id"]: r.get("phenotype_name") for r in runs_wgs if r.get("phenotype_name") in ("Breast Neoplasms", "Health")}

df_ko_long = pd.DataFrame(ko_json)
df_ko_long = df_ko_long[df_ko_long["run_id"].isin(run_ph)]
ko_mat = df_ko_long.pivot_table(index="run_id", columns="ko", values="abundance", aggfunc="sum", fill_value=0.0)
# Filter KOs present in >= 20% samples
ko_mat = ko_mat.loc[:, (ko_mat > 0).mean(axis=0) >= 0.20]
print(f"PRJNA453965 KO matrix: {ko_mat.shape} (BC={sum(run_ph[r]=='Breast Neoplasms' for r in ko_mat.index)}, Health={sum(run_ph[r]=='Health' for r in ko_mat.index)})")

ko_mat.to_csv(os.path.join(hgmt_dir, "HGMT_PRJNA453965_ko_matrix.csv"))

idx_bc_k = np.array([run_ph[r] == "Breast Neoplasms" for r in ko_mat.index])
idx_hc_k = ~idx_bc_k
ko_rows = []
for ko in ko_mat.columns:
    v = ko_mat[ko].values
    v1, v0 = v[idx_bc_k], v[idx_hc_k]
    _, pval = stats.mannwhitneyu(v1, v0, alternative="two-sided")
    ko_rows.append({
        "ko": ko,
        "mean_BC": float(v1.mean()),
        "mean_Health": float(v0.mean()),
        "log2FC": float(np.log2((v1.mean() + 1e-6) / (v0.mean() + 1e-6))),
        "wilcox_p": float(pval),
    })
df_ko_diff = pd.DataFrame(ko_rows)
df_ko_diff["wilcox_fdr"] = bh_fdr(df_ko_diff["wilcox_p"].values)
df_ko_diff = df_ko_diff.sort_values("wilcox_p")
df_ko_diff.to_csv(os.path.join(res_dir, "HGMT_PRJNA453965_WGS_KO_diff.csv"), index=False)
print(f"Significant KOs in PRJNA453965: p<0.05 = {(df_ko_diff['wilcox_p']<0.05).sum()}, FDR<0.05 = {(df_ko_diff['wilcox_fdr']<0.05).sum()}")
print(df_ko_diff.head(20).to_string(index=False))
