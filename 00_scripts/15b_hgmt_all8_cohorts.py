import json
import os
from collections import defaultdict
import numpy as np
import pandas as pd
from scipy import stats

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

hgmt_dir = os.path.join(project_root, "01_data", "HGMT")
res_dir  = os.path.join(project_root, "02_results")

def bh_fdr(pvals):
    p = np.asarray(pvals, dtype=float)
    n = len(p)
    if n == 0:
        return p
    order = np.argsort(p)
    ranks = np.arange(1, n + 1)
    q = np.empty(n)
    q[order] = p[order] * n / ranks
    for i in range(n - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.clip(q, 0, 1)

def parse_genus_matrix(pid, pheno_filter=None):
    runs = json.load(open(os.path.join(hgmt_dir, f"{pid}_runs.json"), encoding="utf-8"))
    bact = json.load(open(os.path.join(hgmt_dir, f"{pid}_feature_bacteria.json"), encoding="utf-8"))
    det  = json.load(open(os.path.join(hgmt_dir, f"{pid}_details.json"), encoding="utf-8"))
    meta_df = pd.DataFrame(runs)
    if pheno_filter is not None:
        meta_df = meta_df[meta_df["phenotype_name"].isin(pheno_filter)].copy()
    keep_set = set(meta_df["run_id"])
    gmap = defaultdict(lambda: defaultdict(float))
    for it in bact:
        rid = it["run_id"]
        if rid not in keep_set:
            continue
        parts = it["taxa"].split("|")
        gp = [x for x in parts if x.startswith("g__")]
        if not gp:
            continue
        gn = gp[-1].replace("g__", "").strip()
        if not gn or gn in ("unclassified", "NA"):
            continue
        gmap[rid][gn] += float(it["abundance"])
    valid_rids = [r for r in meta_df["run_id"] if r in gmap]
    meta_df = meta_df[meta_df["run_id"].isin(valid_rids)].drop_duplicates("run_id").set_index("run_id")
    genera = sorted({g for r in valid_rids for g in gmap[r]})
    mat = pd.DataFrame(0.0, index=valid_rids, columns=genera)
    for r in valid_rids:
        for g, v in gmap[r].items():
            mat.at[r, g] = v
    return det, meta_df, mat

# 1. Summarize ALL 8 HGMT Breast Cancer cohorts
all8 = ["PRJNA345373", "PRJNA383849", "PRJNA453965", "PRJNA494824",
        "PRJNA566060", "PRJNA658160", "PRJNA718520", "PRJNA751260"]
rows_meta = []
all8_data = {}
for pid in all8:
    det, mdf, mat = parse_genus_matrix(pid)
    all8_data[pid] = (det, mdf, mat)
    ph_counts = mdf["phenotype_name"].value_counts().to_dict()
    countries = ",".join(sorted({str(x).strip() for x in mdf["country"].dropna().unique() if str(x) != "NA"}))
    rows_meta.append({
        "project_id": pid,
        "assay_type": det.get("assay_type"),
        "pmid": det.get("public_accession"),
        "country": countries,
        "n_runs_bacteria": len(mat),
        "n_genera": mat.shape[1],
        "n_BC": ph_counts.get("Breast Neoplasms", 0),
        "n_Health": ph_counts.get("Health", 0),
        "n_BreastDiseases": ph_counts.get("Breast Diseases", 0),
        "has_KO": os.path.exists(os.path.join(hgmt_dir, f"{pid}_feature_ko.json")),
        "phenotype_breakdown": "; ".join(f"{k}:{v}" for k, v in ph_counts.items()),
    })
df_meta8 = pd.DataFrame(rows_meta)
df_meta8.to_csv(os.path.join(res_dir, "HGMT_all8_cohorts_metadata_summary.csv"), index=False)
print("=== 1. HGMT All 8 Breast Cancer Cohorts Summary ===")
print(df_meta8.to_string(index=False))

# 2. Stepwise progression (Health -> Breast Diseases -> Breast Neoplasms) in PRJNA658160 & PRJNA751260
print("\n=== 2. Stepwise Dysbiosis (Health -> Benign Breast Diseases -> BC) ===")
grad_rows = []
for pid in ["PRJNA658160", "PRJNA751260"]:
    det, mdf, mat = all8_data[pid]
    mdf = mdf[mdf["phenotype_name"].isin(["Health", "Breast Diseases", "Breast Neoplasms"])].copy()
    mat = mat.loc[mdf.index]
    vals = mat.values
    rs = np.where(vals.sum(axis=1, keepdims=True) > 0, vals.sum(axis=1, keepdims=True), 1.0)
    pr = vals / rs
    with np.errstate(divide="ignore", invalid="ignore"):
        shannon = -np.nansum(np.where(pr > 0, pr * np.log(pr), 0.0), axis=1)
    richness = (vals > 0).sum(axis=1)
    if pid == "PRJNA658160":
        # Local input for Figure 1B, using the same abundance and diversity values.
        plot_input = pd.DataFrame({
            "stage": mdf["phenotype_name"].to_numpy(),
            "Richness": richness,
            "Shannon": shannon,
            "Lachnospira": mat["Lachnospira"].to_numpy(),
            "Roseburia": mat["Roseburia"].to_numpy(),
        })
        plot_input.to_csv(os.path.join(hgmt_dir, f"{pid}_sample_level.csv"), index=False)
    stage_num = mdf["phenotype_name"].map({"Health": 0, "Breast Diseases": 1, "Breast Neoplasms": 2}).values

    for grp in ["Health", "Breast Diseases", "Breast Neoplasms"]:
        idx = (mdf["phenotype_name"] == grp).values
        print(f"  {pid} {grp:17s} n={idx.sum():3d} | Shannon med={np.median(shannon[idx]):.3f} | Richness med={np.median(richness[idx]):.1f}")
    rho_s, p_s = stats.spearmanr(stage_num, shannon)
    rho_r, p_r = stats.spearmanr(stage_num, richness)
    kw_s, p_kw_s = stats.kruskal(shannon[stage_num==0], shannon[stage_num==1], shannon[stage_num==2])
    kw_r, p_kw_r = stats.kruskal(richness[stage_num==0], richness[stage_num==1], richness[stage_num==2])
    print(f"  -> Trend {pid}: Shannon rho={rho_s:.3f} (p={p_s:.3e}, KW p={p_kw_s:.3e}) | Richness rho={rho_r:.3f} (p={p_r:.3e}, KW p={p_kw_r:.3e})")

    for g in ["Lachnospira", "Roseburia", "Ruminococcus_E", "Parasutterella", "Turicibacter", "Holdemanella", "Eggerthella"]:
        if g not in mat.columns:
            continue
        v = mat[g].values
        rho_g, p_g = stats.spearmanr(stage_num, v)
        grad_rows.append({
            "project_id": pid, "feature": g, "spearman_rho": rho_g, "spearman_p": p_g,
            "mean_Health": v[stage_num==0].mean(), "mean_BreastDiseases": v[stage_num==1].mean(), "mean_BC": v[stage_num==2].mean()
        })
    grad_rows.append({
        "project_id": pid, "feature": "Shannon", "spearman_rho": rho_s, "spearman_p": p_s,
        "mean_Health": shannon[stage_num==0].mean(), "mean_BreastDiseases": shannon[stage_num==1].mean(), "mean_BC": shannon[stage_num==2].mean()
    })
    grad_rows.append({
        "project_id": pid, "feature": "Richness", "spearman_rho": rho_r, "spearman_p": p_r,
        "mean_Health": richness[stage_num==0].mean(), "mean_BreastDiseases": richness[stage_num==1].mean(), "mean_BC": richness[stage_num==2].mean()
    })
pd.DataFrame(grad_rows).to_csv(os.path.join(res_dir, "HGMT_benign_to_malignant_gradient.csv"), index=False)

# 3. Menopause-stratified WGS analysis in PRJNA453965 (Pre vs Post) + external WGS comparison (PRJNA718520, PRJNA494824)
print("\n=== 3. PRJNA453965 WGS Menopause-Stratified Analysis ===")
det453, mdf453, mat453 = all8_data["PRJNA453965"]
ko453 = pd.read_csv(os.path.join(hgmt_dir, "HGMT_PRJNA453965_ko_matrix.csv"), index_col=0)
pway453 = pd.read_csv(os.path.join(hgmt_dir, "HGMT_PRJNA453965_pathway_matrix.csv"), index_col=0)

print("PRJNA453965 phenotype x menopause:")
print(pd.crosstab(mdf453["phenotype_name"], mdf453["sample_description"]))

meno_rows = []
for stratum in ["Postmenopause", "Premenopause"]:
    sub_idx = mdf453.index[mdf453["sample_description"] == stratum]
    idx_bc = [r for r in sub_idx if mdf453.loc[r, "phenotype_name"] == "Breast Neoplasms"]
    idx_hc = [r for r in sub_idx if mdf453.loc[r, "phenotype_name"] == "Health"]
    for g in ["Lachnospira", "Roseburia", "Escherichia", "Turicibacter", "Lawsonibacter", "Eggerthella"]:
        if g not in mat453.columns:
            continue
        v1 = mat453.loc[idx_bc, g].values
        v0 = mat453.loc[idx_hc, g].values
        _, pv = stats.mannwhitneyu(v1, v0, alternative="two-sided")
        meno_rows.append({
            "stratum": stratum, "type": "Genus", "feature": g, "n_BC": len(idx_bc), "n_Health": len(idx_hc),
            "mean_BC": v1.mean(), "mean_Health": v0.mean(),
            "log2FC": np.log2((v1.mean()+1e-4)/(v0.mean()+1e-4)), "wilcox_p": pv
        })
    for pw in ["map00984", "map00380", "map00071", "map00650", "map00220", "map01040"]:
        if pw not in pway453.columns:
            continue
        v1 = pway453.loc[[r for r in idx_bc if r in pway453.index], pw].values
        v0 = pway453.loc[[r for r in idx_hc if r in pway453.index], pw].values
        _, pv = stats.mannwhitneyu(v1, v0, alternative="two-sided")
        meno_rows.append({
            "stratum": stratum, "type": "KEGG_Pathway", "feature": pw, "n_BC": len(v1), "n_Health": len(v0),
            "mean_BC": v1.mean(), "mean_Health": v0.mean(),
            "log2FC": np.log2((v1.mean()+1e-6)/(v0.mean()+1e-6)), "wilcox_p": pv
        })
df_meno = pd.DataFrame(meno_rows)
df_meno.to_csv(os.path.join(res_dir, "HGMT_PRJNA453965_WGS_menopause_stratified.csv"), index=False)
print(df_meno.to_string(index=False))

# 4. External WGS BC cohorts (PRJNA718520 France n=76 BC, PRJNA494824 n=8 BC) pathway & genus abundance
print("\n=== 4. External WGS Breast Cancer Cohorts (PRJNA718520 & PRJNA494824) ===")
ko_link = pd.read_csv(os.path.join(hgmt_dir, "kegg_ko_pathway_link.tsv"), sep="\t", names=["ko_raw", "pw_raw"])
ko_link["ko"] = ko_link["ko_raw"].str.replace("ko:", "", regex=False)
ko_link["pathway"] = ko_link["pw_raw"].str.replace("path:", "", regex=False)
ko_link = ko_link[ko_link["pathway"].str.startswith("map")]
pw_names = pd.read_csv(os.path.join(hgmt_dir, "kegg_pathway_names.tsv"), sep="\t", names=["pw_raw", "name"])
pw_names["pathway"] = pw_names["pw_raw"].str.replace("path:", "", regex=False)
pw_map = dict(zip(pw_names["pathway"], pw_names["name"]))

wgs_comp_rows = []
for pid in ["PRJNA453965", "PRJNA718520", "PRJNA494824"]:
    det, mdf, mat = all8_data[pid]
    bc_rids = mdf.index[mdf["phenotype_name"] == "Breast Neoplasms"].tolist()
    # Save genus table for non-CC cohorts too
    if pid in ("PRJNA718520", "PRJNA494824", "PRJNA566060"):
        mout = mat.copy()
        mout.insert(0, "phenotype", mdf.loc[mat.index, "phenotype_name"])
        mout.insert(0, "project_id", pid)
        mout.to_csv(os.path.join(hgmt_dir, f"HGMT_genus_{pid}.csv"))

    # Compute KO pathway abundances for BC samples
    if pid == "PRJNA453965":
        pw_bc = pway453.loc[[r for r in bc_rids if r in pway453.index]]
        pw_hc = pway453.loc[mdf.index[mdf["phenotype_name"] == "Health"].tolist()]
    else:
        kj = json.load(open(os.path.join(hgmt_dir, f"{pid}_feature_ko.json"), encoding="utf-8"))
        dk = pd.DataFrame(kj)
        dk = dk[dk["run_id"].isin(bc_rids)]
        km = dk.pivot_table(index="run_id", columns="ko", values="abundance", aggfunc="sum", fill_value=0.0)
        sub_link = ko_link[ko_link["ko"].isin(km.columns)]
        pw_dict = {}
        for pw, grp in sub_link.groupby("pathway"):
            cols = [c for c in grp["ko"].unique() if c in km.columns]
            if len(cols) >= 3:
                pw_dict[pw] = km[cols].sum(axis=1)
        pw_bc = pd.DataFrame(pw_dict)
        pw_bc.to_csv(os.path.join(hgmt_dir, f"HGMT_{pid}_BC_pathway_matrix.csv"))

    for g in ["Lachnospira", "Roseburia", "Escherichia", "Faecalibacterium", "Turicibacter", "Bacteroides", "Eggerthella"]:
        v = mat.loc[bc_rids, g].values if g in mat.columns else np.zeros(len(bc_rids))
        wgs_comp_rows.append({
            "project_id": pid, "group": "BC", "n_samples": len(bc_rids),
            "feature_type": "Genus", "feature": g, "mean": v.mean(), "prevalence": (v > 0).mean()
        })
    for pw in ["map00984", "map00380", "map00071", "map00650", "map00220", "map01040"]:
        v = pw_bc[pw].values if pw in pw_bc.columns else np.zeros(len(bc_rids))
        wgs_comp_rows.append({
            "project_id": pid, "group": "BC", "n_samples": len(pw_bc),
            "feature_type": "KEGG_Pathway", "feature": pw, "mean": v.mean(), "prevalence": (v > 0).mean()
        })
df_wcomp = pd.DataFrame(wgs_comp_rows)
df_wcomp.to_csv(os.path.join(res_dir, "HGMT_3WGS_BC_cross_cohort_profile.csv"), index=False)
print("Saved HGMT_3WGS_BC_cross_cohort_profile.csv; top summary:")
print(df_wcomp.to_string(index=False))
