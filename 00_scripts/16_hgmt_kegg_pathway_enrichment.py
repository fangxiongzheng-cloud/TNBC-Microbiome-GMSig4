import os
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
    q = np.empty(n)
    q[order] = p[order] * n / np.arange(1, n + 1)
    for i in range(n - 2, -1, -1):
        q[order[i]] = min(q[order[i]], q[order[i + 1]])
    return np.clip(q, 0, 1)

# Load KEGG mappings
ko_pw = pd.read_csv(os.path.join(hgmt_dir, "kegg_ko_pathway_link.tsv"), sep="\t", header=None, names=["ko", "pathway"])
ko_pw["ko"] = ko_pw["ko"].str.replace("ko:", "", regex=False)
ko_pw["pathway"] = ko_pw["pathway"].str.replace("path:", "", regex=False)
# Keep map pathways (avoid duplicate ko/map pairs: keep mapXXXXX)
ko_pw = ko_pw[ko_pw["pathway"].str.startswith("map")].drop_duplicates()

pw_names = pd.read_csv(os.path.join(hgmt_dir, "kegg_pathway_names.tsv"), sep="\t", header=None, names=["pathway", "pathway_name"])
pw_names["pathway"] = pw_names["pathway"].str.replace("path:", "", regex=False).str.replace("^ko", "map", regex=True)

ko_names = pd.read_csv(os.path.join(hgmt_dir, "kegg_ko_names.tsv"), sep="\t", header=None, names=["ko", "ko_name"])
ko_names["ko"] = ko_names["ko"].str.replace("ko:", "", regex=False)

# Annotate PRJNA453965 KO diff table
ko_diff = pd.read_csv(os.path.join(res_dir, "HGMT_PRJNA453965_WGS_KO_diff.csv"))
ko_diff_ann = ko_diff.merge(ko_names, on="ko", how="left")
ko_diff_ann.to_csv(os.path.join(res_dir, "HGMT_PRJNA453965_WGS_KO_diff_annotated.csv"), index=False)

# Pathway-level abundance aggregation across all 133 WGS samples (62 BC vs 71 Health)
import json
runs_wgs = json.load(open(os.path.join(hgmt_dir, "PRJNA453965_runs.json"), encoding="utf-8"))
run_ph = {r["run_id"]: r.get("phenotype_name") for r in runs_wgs if r.get("phenotype_name") in ("Breast Neoplasms", "Health")}

ko_mat = pd.read_csv(os.path.join(hgmt_dir, "HGMT_PRJNA453965_ko_matrix.csv"), index_col=0)
# Map KOs in ko_mat to pathways
sub_map = ko_pw[ko_pw["ko"].isin(ko_mat.columns)]
pw_abund = {}
for pw, grp in sub_map.groupby("pathway"):
    kos = grp["ko"].unique()
    if len(kos) < 3:
        continue
    pw_abund[pw] = ko_mat[kos].sum(axis=1)

df_pw_mat = pd.DataFrame(pw_abund, index=ko_mat.index)
df_pw_mat.to_csv(os.path.join(hgmt_dir, "HGMT_PRJNA453965_pathway_matrix.csv"))

idx_bc = np.array([run_ph[r] == "Breast Neoplasms" for r in df_pw_mat.index])
idx_hc = ~idx_bc

rows = []
for pw in df_pw_mat.columns:
    v = df_pw_mat[pw].values
    v1, v0 = v[idx_bc], v[idx_hc]
    _, pval = stats.mannwhitneyu(v1, v0, alternative="two-sided")
    rows.append({
        "pathway": pw,
        "n_kos": int((sub_map["pathway"] == pw).sum()),
        "mean_BC": float(v1.mean()),
        "mean_Health": float(v0.mean()),
        "log2FC": float(np.log2((v1.mean() + 1e-6) / (v0.mean() + 1e-6))),
        "wilcox_p": float(pval),
    })

df_pw_diff = pd.DataFrame(rows).merge(pw_names, on="pathway", how="left")
df_pw_diff["wilcox_fdr"] = bh_fdr(df_pw_diff["wilcox_p"].values)
df_pw_diff = df_pw_diff.sort_values("wilcox_p")
df_pw_diff.to_csv(os.path.join(res_dir, "HGMT_PRJNA453965_WGS_KEGG_pathway_diff.csv"), index=False)

print(f"Total pathways tested (>=3 KOs): {len(df_pw_diff)}, FDR<0.05: {(df_pw_diff['wilcox_fdr']<0.05).sum()}")
print("\nTop 20 differential KEGG pathways in WGS Breast Cancer:")
print(df_pw_diff.head(20)[["pathway", "pathway_name", "n_kos", "log2FC", "wilcox_p", "wilcox_fdr"]].to_string(index=False))

# Check key microbiome-immune pathways explicitly
key_terms = ["butanoate", "propanoate", "tryptophan", "bile", "lipopolysaccharide", "fatty acid", "phenylalanine", "tyrosine", "steroid", "estrogen", "folate"]
mask = df_pw_diff["pathway_name"].fillna("").str.lower().apply(lambda x: any(t in x for t in key_terms))
print("\n=== Microbiome-Metabolic Relevant KEGG Pathways ===")
print(df_pw_diff[mask][["pathway", "pathway_name", "n_kos", "log2FC", "wilcox_p", "wilcox_fdr"]].to_string(index=False))
