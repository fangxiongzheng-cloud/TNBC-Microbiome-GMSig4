import os
import pandas as pd
import numpy as np
import scanpy as sc
from scipy.stats import spearmanr

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

out_h5ad = os.path.join(project_root, "01_data", "scRNA", "GSE176078_TNBC_annotated.h5ad")
res_dir  = os.path.join(project_root, "02_results")

adata = sc.read_h5ad(out_h5ad)
print(f"Loaded cached TNBC atlas: {adata.n_obs} cells x {adata.n_vars} genes", flush=True)

gut_df = pd.read_csv(os.path.join(res_dir, "gutMGene2_human_target_genes_master.csv"))
gut_genes = [g for g in gut_df["gene"].unique() if g in adata.var_names]

X_gut = adata[:, gut_genes].X
if hasattr(X_gut, "toarray"):
    X_gut = X_gut.toarray()
expr_df = pd.DataFrame(X_gut, index=adata.obs_names, columns=gut_genes)

pb = expr_df.groupby(adata.obs["orig.ident"].astype(str), observed=True)[gut_genes].mean()
ct_prop = pd.crosstab(adata.obs["orig.ident"].astype(str), adata.obs["celltype_major"].astype(str), normalize="index")
ct_prop.columns = [f"prop_{c}" for c in ct_prop.columns]
pb_full = pd.concat([pb, ct_prop], axis=1)
pb_full.to_csv(os.path.join(res_dir, "scRNA_TNBC_patient_pseudobulk_and_props.csv"))

# Also minor subtype T-cell / Myeloid specificity for signature + causal driver genes
mean_major = pd.read_csv(os.path.join(res_dir, "scRNA_TNBC_gutMGene_celltype_major_mean.csv"), index_col=0)
pct_major  = pd.read_csv(os.path.join(res_dir, "scRNA_TNBC_gutMGene_celltype_major_pct.csv"), index_col=0)
mean_minor = pd.read_csv(os.path.join(res_dir, "scRNA_TNBC_gutMGene_celltype_minor_mean.csv"), index_col=0)

focus_genes = [g for g in ["IFNG", "IDO1", "CXCL10", "ENO1", "TGFB1", "TGFBR1", "TCF7L2", "PEMT",
                           "MAGI3", "CCR6", "FOXP3", "CXCR6", "LYZ", "CASP1", "CLDN7", "F11R",
                           "ST6GAL1", "LDHA", "DAPK3", "VEGFA", "MUC1"] if g in mean_major.index]
print("\n=== Single-cell Major Cell-type Mean Expression (TNBC n=42,512 cells, 10 patients) ===")
print(mean_major.loc[focus_genes].round(3).to_string())

print("\n=== Single-cell % Expressing Cells per Major Cell Type ===")
print((pct_major.loc[focus_genes] * 100).round(1).to_string())

# Also compute C16 and C1 signature pseudobulk risk score in the 10 TNBC patients (z-scored across 10 patients)
# Using exact TCGA-TNBC UniBeta coefficients from Signature_Candidate_Selection_All5Cohorts.csv
c16_coef = {"IFNG": -0.5941, "IDO1": -0.4378, "CXCL10": -0.3292, "ENO1": 0.2195}
c1_coef  = {"IFNG": -0.5941, "IDO1": -0.4378, "CXCL10": -0.3292, "TGFB1": 0.6052}

pb_z = (pb[list(set(c16_coef) | set(c1_coef))] - pb[list(set(c16_coef) | set(c1_coef))].mean()) / pb[list(set(c16_coef) | set(c1_coef))].std(ddof=1)
pb_full["Risk_C16_IfngIdo1Cxcl10Eno1"] = sum(pb_z[g] * w for g, w in c16_coef.items())
pb_full["Risk_C1_IfngIdo1Cxcl10Tgfb1"] = sum(pb_z[g] * w for g, w in c1_coef.items())

# T-cell CD8+ Exhausted / Cytotoxic minor subtype proportions per patient
minor_prop = pd.crosstab(adata.obs["orig.ident"].astype(str), adata.obs["celltype_minor"].astype(str), normalize="index")
pb_full = pd.concat([pb_full, minor_prop.add_prefix("minor_prop_")], axis=1)
pb_full.to_csv(os.path.join(res_dir, "scRNA_TNBC_patient_pseudobulk_and_props.csv"))

print("\n=== Correlation between Patient Pseudobulk Risk Score (n=10 TNBC tumors) and TME Cell Proportions ===")
for rcol in ["Risk_C16_IfngIdo1Cxcl10Eno1", "Risk_C1_IfngIdo1Cxcl10Tgfb1"]:
    print(f"\n-- {rcol} --")
    for pcol in [c for c in pb_full.columns if c.startswith("prop_") or c.startswith("minor_prop_T") or c.startswith("minor_prop_Cancer")]:
        rho, pv = spearmanr(pb_full[rcol], pb_full[pcol])
        if pv < 0.20:
            print(f"  {pcol:35s} rho={rho:+.3f}, p={pv:.4g}")

# Also compute C16 & C1 single-cell module scores inside T-cells vs Myeloid vs Epithelial vs CAFs
cell_focus = mean_minor.loc[focus_genes, [c for c in mean_minor.columns if any(k in c for k in ["T_cells", "Myeloid", "Cancer", "CAF", "Plasmablasts", "B_cells"])]]
cell_focus.to_csv(os.path.join(res_dir, "scRNA_TNBC_focus_genes_minor_celltypes.csv"))
print("\nDONE scRNA TME correlation & focus gene export.", flush=True)
