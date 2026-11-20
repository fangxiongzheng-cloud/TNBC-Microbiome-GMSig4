import os
import pandas as pd
import numpy as np
import scipy.io
import scanpy as sc

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")

sc.settings.verbosity = 2

scrna_dir = os.path.join(project_root, "01_data", "scRNA", "Wu_etal_2021_BRCA_scRNASeq")
out_h5ad  = os.path.join(project_root, "01_data", "scRNA", "GSE176078_TNBC_annotated.h5ad")
res_dir   = os.path.join(project_root, "02_results")

print("1. Loading metadata, genes, barcodes...", flush=True)
meta = pd.read_csv(os.path.join(scrna_dir, "metadata.csv"), index_col=0)
genes = pd.read_csv(os.path.join(scrna_dir, "count_matrix_genes.tsv"), header=None, names=["gene"])["gene"].astype(str).values
barcodes = pd.read_csv(os.path.join(scrna_dir, "count_matrix_barcodes.tsv"), header=None, names=["barcode"])["barcode"].astype(str).values

print("Subtype breakdown in full BRCA atlas (n=100,064):", flush=True)
print(meta["subtype"].value_counts(), flush=True)

tnbc_mask = (meta.loc[barcodes, "subtype"].values == "TNBC")
tnbc_idx = np.where(tnbc_mask)[0]
print(f"TNBC cells: {len(tnbc_idx)} across {meta.loc[barcodes[tnbc_mask], 'orig.ident'].nunique()} patients", flush=True)
print(meta.loc[barcodes[tnbc_mask], "orig.ident"].value_counts(), flush=True)
print(meta.loc[barcodes[tnbc_mask], "celltype_major"].value_counts(), flush=True)

print("2. Reading sparse count matrix (2.3 GB) in CSC/CSR...", flush=True)
mtx = scipy.io.mmread(os.path.join(scrna_dir, "count_matrix_sparse.mtx")).tocsc()
print(f"Loaded matrix shape: {mtx.shape}, slicing TNBC columns...", flush=True)
mtx_tnbc = mtx[:, tnbc_idx].T.tocsr()  # cells x genes
del mtx

tnbc_barcodes = barcodes[tnbc_idx]
tnbc_meta = meta.loc[tnbc_barcodes].copy()

adata = sc.AnnData(X=mtx_tnbc, obs=tnbc_meta, var=pd.DataFrame(index=genes))
adata.var_names_make_unique()
adata.layers["counts"] = adata.X.copy()

print("3. Normalizing (CP10K + log1p) and running HVG/PCA/UMAP on TNBC atlas...", flush=True)
sc.pp.normalize_total(adata, target_sum=1e4)
sc.pp.log1p(adata)
adata.raw = adata

sc.pp.highly_variable_genes(adata, min_mean=0.0125, max_mean=3, min_disp=0.5, batch_key="orig.ident")
adata_hvg = adata[:, adata.var["highly_variable"]].copy()
sc.pp.scale(adata_hvg, max_value=10)
sc.tl.pca(adata_hvg, n_comps=30, svd_solver="arpack")
sc.pp.neighbors(adata_hvg, n_neighbors=15, n_pcs=30)
sc.tl.umap(adata_hvg, random_state=42)

adata.obsm["X_pca"] = adata_hvg.obsm["X_pca"]
adata.obsm["X_umap"] = adata_hvg.obsm["X_umap"]

adata.write_h5ad(out_h5ad, compression="gzip")
print(f"SAVED TNBC AnnData to {out_h5ad} ({adata.n_obs} cells x {adata.n_vars} genes)", flush=True)

# 4. Compute celltype composition per TNBC patient and gutMGene cell-type specificity
gut_df = pd.read_csv(os.path.join(res_dir, "gutMGene2_human_target_genes_master.csv"))
gut_genes = [g for g in gut_df["gene"].unique() if g in adata.var_names]
print(f"gutMGene genes detected in scRNA-seq: {len(gut_genes)} / {gut_df['gene'].nunique()}", flush=True)

X_gut = adata[:, gut_genes].X
if hasattr(X_gut, "toarray"):
    X_gut = X_gut.toarray()
expr_df = pd.DataFrame(X_gut, index=adata.obs_names, columns=gut_genes)
expr_df["celltype_major"] = adata.obs["celltype_major"].values
expr_df["celltype_minor"] = adata.obs["celltype_minor"].values
expr_df["patient"] = adata.obs["orig.ident"].values

mean_major = expr_df.groupby("celltype_major")[gut_genes].mean()
mean_major.T.to_csv(os.path.join(res_dir, "scRNA_TNBC_gutMGene_celltype_major_mean.csv"))

pct_major = expr_df.groupby("celltype_major")[gut_genes].apply(lambda d: (d > 0).mean())
pct_major.T.to_csv(os.path.join(res_dir, "scRNA_TNBC_gutMGene_celltype_major_pct.csv"))

mean_minor = expr_df.groupby("celltype_minor")[gut_genes].mean()
mean_minor.T.to_csv(os.path.join(res_dir, "scRNA_TNBC_gutMGene_celltype_minor_mean.csv"))

# Save UMAP coords + major/minor/subset annotations + key genes for clean Figure plotting
key_cols = [g for g in ["IFNG", "IDO1", "CXCL10", "ENO1", "TGFB1", "TGFBR1", "MAGI3", "TCF7L2", "CCR6", "CXCR6", "LYZ", "CASP1", "CLDN7", "F11R", "VEGFA", "DAPK3", "CD8A", "CD274", "PDCD1", "FOXP3"] if g in expr_df.columns]
umap_df = pd.DataFrame({
    "cell_id": adata.obs_names,
    "patient": adata.obs["orig.ident"].values,
    "celltype_major": adata.obs["celltype_major"].values,
    "celltype_minor": adata.obs["celltype_minor"].values,
    "celltype_subset": adata.obs["celltype_subset"].values,
    "UMAP_1": adata.obsm["X_umap"][:, 0],
    "UMAP_2": adata.obsm["X_umap"][:, 1],
})
for g in key_cols:
    umap_df[g] = expr_df[g].values
umap_df.to_csv(os.path.join(res_dir, "scRNA_TNBC_UMAP_metadata_keygenes.csv"), index=False)

# Patient-level pseudobulk expression + TME proportions
pb = expr_df.groupby("patient")[gut_genes].mean()
ct_prop = pd.crosstab(adata.obs["orig.ident"], adata.obs["celltype_major"], normalize="index")
pb.join(ct_prop, lsuffix="", rsuffix="_prop").to_csv(os.path.join(res_dir, "scRNA_TNBC_patient_pseudobulk_and_props.csv"))
print("DONE scRNA TNBC atlas extraction & cell-type specificity tables!", flush=True)
