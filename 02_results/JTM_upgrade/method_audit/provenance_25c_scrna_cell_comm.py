import os
import re
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import stats
import matplotlib

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec

DATA_DIR = os.path.join(project_root, "01_data")
RES_DIR  = os.path.join(project_root, "02_results")
SUP_DIR  = os.path.join(project_root, "03_figures", "Supplementary")
os.makedirs(SUP_DIR, exist_ok=True)

print("=== 1. Loading Annotated 42,512-Cell TNBC Atlas ===", flush=True)
adata = sc.read_h5ad(os.path.join(DATA_DIR, "scRNA", "GSE176078_TNBC_annotated.h5ad"))
pb_df = pd.read_csv(os.path.join(RES_DIR, "scRNA_TNBC_patient_pseudobulk_and_props.csv"))

# Stratify 10 TNBC patients into High vs Low GMSig-4 Risk (5 vs 5 median split)
pb_df = pb_df.sort_values("Risk_C16_IfngIdo1Cxcl10Eno1")
med_risk = pb_df["Risk_C16_IfngIdo1Cxcl10Eno1"].median()
pb_df["sc_RiskGroup"] = np.where(pb_df["Risk_C16_IfngIdo1Cxcl10Eno1"] > med_risk, "High_GMSig4", "Low_GMSig4")
pt_grp_map = dict(zip(pb_df["orig.ident"], pb_df["sc_RiskGroup"]))
print("Patient GMSig-4 single-cell stratification (5 High vs 5 Low):", flush=True)
print(pb_df[["orig.ident", "Risk_C16_IfngIdo1Cxcl10Eno1", "prop_T-cells", "prop_Myeloid", "prop_Cancer Epithelial", "sc_RiskGroup"]].to_string(index=False), flush=True)

adata.obs["sc_RiskGroup"] = adata.obs["orig.ident"].map(pt_grp_map).astype(str)
# Clean cell niche grouping (merging Plasmablasts+B-cells as B/Plasma only if tiny, or keep clear functional niches; keep all main immune/stromal/tumor!)
# Check minor subsets for CD4/CD8/Macrophage/DC/iCAF/myCAF/Basal/Cycling richness
niche_col = "celltype_major"
niches = ["Cancer Epithelial", "CAFs", "Endothelial", "Myeloid", "T-cells", "B-cells", "Plasmablasts", "PVL", "Normal Epithelial"]
# Also build focused high-resolution TME sub-niches for detailed chemokine/checkpoint axis
sub_map = {
    "T cells CD8+": "CD8+ T",
    "T cells CD4+": "CD4+ T",
    "Cycling T-cells": "Cycling T",
    "NK cells": "NK/NKT",
    "NKT cells": "NK/NKT",
    "Macrophage": "Macrophage",
    "Monocyte": "Monocyte",
    "DCs": "DCs",
    "Cycling_Myeloid": "Cycling Myeloid",
    "CAFs MSC iCAF-like": "iCAF",
    "CAFs myCAF-like": "myCAF",
    "Cancer Basal SC": "Cancer Basal",
    "Cancer Cycling": "Cancer Cycling",
    "Endothelial ACKR1": "Endothelial",
    "Endothelial CXCL12": "Endothelial",
    "Endothelial RGS5": "Endothelial",
    "Endothelial Lymphatic LYVE1": "Endothelial",
}
adata.obs["TME_niche"] = adata.obs["celltype_minor"].astype(str).map(sub_map).fillna(adata.obs["celltype_major"].astype(str))

# Curated Immune Checkpoint, Chemokine Recruitment, Metabolic & TGFB/Microbiome-Linked LR Database
LR_PAIRS = [
    ("CXCL10",  "CXCR3",  "Myeloid-T Axis"),
    ("CXCL9",   "CXCR3",  "Myeloid-T Axis"),
    ("CXCL11",  "CXCR3",  "Myeloid-T Axis"),
    ("CXCL16",  "CXCR6",  "Myeloid-T Axis"),
    ("CCL5",    "CCR5",   "T/Myeloid Chemokine"),
    ("CCL4",    "CCR5",   "T/Myeloid Chemokine"),
    ("CCL3",    "CCR5",   "T/Myeloid Chemokine"),
    ("CCL20",   "CCR6",   "Microbiome-B/T Axis"),
    ("CXCL12",  "CXCR4",  "CAF-Immune Axis"),
    ("CCL2",    "CCR2",   "CAF/Myeloid Recruit"),
    ("IFNG",    "IFNGR1", "Effector IFN-gamma"),
    ("IFNG",    "IFNGR2", "Effector IFN-gamma"),
    ("CD274",   "PDCD1",  "Immune Checkpoint"),
    ("PDCD1LG2","PDCD1",  "Immune Checkpoint"),
    ("LGALS9",  "HAVCR2", "Immune Checkpoint"),
    ("NECTIN2", "TIGIT",  "Immune Checkpoint"),
    ("PVR",     "TIGIT",  "Immune Checkpoint"),
    ("CD86",    "CTLA4",  "Immune Checkpoint"),
    ("CD80",    "CD28",   "Co-stimulation"),
    ("CD86",    "CD28",   "Co-stimulation"),
    ("TGFB1",   "TGFBR1", "Stromal Exclusion"),
    ("TGFB1",   "TGFBR2", "Stromal Exclusion"),
    ("MIF",     "CD74",   "Tumor-Myeloid MIF"),
    ("MIF",     "CXCR4",  "Tumor-Immune MIF"),
    ("SPP1",    "CD44",   "Macrophage-CAF/Tumor"),
    ("VEGFA",   "FLT1",   "Angiogenesis"),
    ("VEGFA",   "KDR",    "Angiogenesis"),
    ("COL1A1",  "CD44",   "CAF ECM Exclusion"),
    ("FN1",     "CD44",   "CAF ECM Exclusion"),
    ("FN1",     "ITGB1",  "CAF ECM Exclusion"),
    ("TNF",     "TNFRSF1A","Inflammatory Cytokine"),
    ("ICAM1",   "ITGAL",  "Adhesion / T-Priming"),
]

all_genes = sorted(list({g for l, r, _ in LR_PAIRS for g in (l, r) if g in adata.var_names}))
X_sub = adata[:, all_genes].X
if hasattr(X_sub, "toarray"):
    X_sub = X_sub.toarray()
df_ex = pd.DataFrame(X_sub, index=adata.obs_names, columns=all_genes)
df_ex["celltype_major"] = adata.obs["celltype_major"].astype(str).values
df_ex["TME_niche"]      = adata.obs["TME_niche"].astype(str).values
df_ex["sc_RiskGroup"]   = adata.obs["sc_RiskGroup"].astype(str).values
df_ex["patient"]        = adata.obs["orig.ident"].astype(str).values

# --- Part A: Major Cell-Type Empirical Permutation Cell-Cell Communication (1,000 perms, overall + Low vs High GMSig-4) ---
def compute_comm(df_in, ct_col, lr_list, n_perm=1000, seed=42):
    rng = np.random.default_rng(seed)
    cts = sorted(df_in[ct_col].unique())
    g_mean = df_in.groupby(ct_col)[all_genes].mean()
    g_pct  = df_in.groupby(ct_col)[all_genes].apply(lambda d: (d > 0.1).mean())

    ct_labels = df_in[ct_col].values
    n_cells = len(ct_labels)
    ct_idx_map = {ct: np.where(ct_labels == ct)[0] for ct in cts}
    # Pre-build permutation mean null distribution per gene across cell types!
    X_mat = df_in[all_genes].values
    perm_means = {ct: np.zeros((n_perm, len(all_genes)), dtype=np.float32) for ct in cts}
    sizes = [len(ct_idx_map[ct]) for ct in cts]
    for p in range(n_perm):
        shuf = rng.permutation(n_cells)
        pos = 0
        for ct, sz in zip(cts, sizes):
            perm_means[ct][p, :] = X_mat[shuf[pos:pos+sz]].mean(axis=0)
            pos += sz
    gene_idx = {g: i for i, g in enumerate(all_genes)}

    rows = []
    for lig, rec, axis in lr_list:
        if lig not in gene_idx or rec not in gene_idx:
            continue
        li, ri = gene_idx[lig], gene_idx[rec]
        for s in cts:
            ml = float(g_mean.loc[s, lig]); pl = float(g_pct.loc[s, lig])
            if ml <= 0.02 or pl < 0.03:
                continue
            for t in cts:
                mr = float(g_mean.loc[t, rec]); pr = float(g_pct.loc[t, rec])
                if mr <= 0.02 or pr < 0.03:
                    continue
                score = float(np.sqrt(ml * mr))
                null_scores = np.sqrt(perm_means[s][:, li] * perm_means[t][:, ri])
                pval = float((np.sum(null_scores >= score) + 1.0) / (n_perm + 1.0))
                rows.append((s, t, lig, rec, f"{lig}->{rec}", axis, ml, mr, pl, pr, score, pval))
    out = pd.DataFrame(rows, columns=["source", "target", "ligand", "receptor", "lr_pair", "pathway",
                                      "mean_lig", "mean_rec", "pct_lig", "pct_rec", "comm_score", "perm_p"])
    if len(out) > 0:
        # Benjamini-Hochberg FDR
        rk = out["perm_p"].rank(method="min")
        out["perm_fdr"] = np.minimum(1.0, out["perm_p"] * len(out) / rk)
    return out

print("2. Computing Overall Major-Celltype Permutation Cell Communication (1,000 perms)...", flush=True)
comm_all = compute_comm(df_ex, "celltype_major", LR_PAIRS, n_perm=1000, seed=42)
comm_all.to_csv(os.path.join(RES_DIR, "scRNA_TNBC_CellComm_Interactions.csv"), index=False)
print(f"Total active major interactions tested: {len(comm_all)}, FDR<0.05: {(comm_all['perm_fdr']<0.05).sum()}", flush=True)

print("3. Comparing Low-GMSig4 (n=5) vs High-GMSig4 (n=5) Patients across Major & TME Sub-Niches...", flush=True)
comm_low  = compute_comm(df_ex[df_ex["sc_RiskGroup"] == "Low_GMSig4"],  "celltype_major", LR_PAIRS, n_perm=500, seed=101)
comm_high = compute_comm(df_ex[df_ex["sc_RiskGroup"] == "High_GMSig4"], "celltype_major", LR_PAIRS, n_perm=500, seed=202)
comp = comm_all[["source", "target", "lr_pair", "pathway", "comm_score", "perm_fdr"]].merge(
    comm_low[["source", "target", "lr_pair", "comm_score", "perm_p"]].rename(columns={"comm_score": "score_Low", "perm_p": "p_Low"}),
    on=["source", "target", "lr_pair"], how="left"
).merge(
    comm_high[["source", "target", "lr_pair", "comm_score", "perm_p"]].rename(columns={"comm_score": "score_High", "perm_p": "p_High"}),
    on=["source", "target", "lr_pair"], how="left"
).fillna({"score_Low": 0.0, "score_High": 0.0, "p_Low": 1.0, "p_High": 1.0})

comp["diff_HighMinusLow"] = comp["score_High"] - comp["score_Low"]
comp["log2FC_HighVsLow"]  = np.log2((comp["score_High"] + 0.01) / (comp["score_Low"] + 0.01))
comp.to_csv(os.path.join(RES_DIR, "scRNA_TNBC_CellComm_HighVsLow_GMSig4.csv"), index=False)

# Patient-level interaction score per key LR pair across 10 TNBC patients vs GMSig-4 Risk Score
pt_rows = []
for pt, psub in df_ex.groupby("patient"):
    risk_val = float(pb_df.loc[pb_df["orig.ident"] == pt, "Risk_C16_IfngIdo1Cxcl10Eno1"].iloc[0])
    rgrp     = pt_grp_map[pt]
    gm = psub.groupby("celltype_major")[all_genes].mean()
    for lig, rec, s_ct, t_ct in [
        ("CXCL10", "CXCR3", "Myeloid", "T-cells"),
        ("CXCL9",  "CXCR3", "Myeloid", "T-cells"),
        ("IFNG",   "IFNGR1","T-cells", "Myeloid"),
        ("IFNG",   "IFNGR2","T-cells", "Myeloid"),
        ("CD274",  "PDCD1", "Myeloid", "T-cells"),
        ("CCL5",   "CCR5",  "T-cells", "Myeloid"),
        ("CCL20",  "CCR6",  "Myeloid", "B-cells"),
        ("CXCL16", "CXCR6", "Myeloid", "T-cells"),
        ("MIF",    "CD74",  "Cancer Epithelial", "Myeloid"),
        ("SPP1",   "CD44",  "Myeloid", "Cancer Epithelial"),
        ("COL1A1", "CD44",  "CAFs",    "T-cells"),
        ("TGFB1",  "TGFBR1","Myeloid", "CAFs"),
        ("VEGFA",  "KDR",   "Cancer Epithelial", "Endothelial"),
    ]:
        if s_ct in gm.index and t_ct in gm.index:
            sc_v = float(np.sqrt(gm.loc[s_ct, lig] * gm.loc[t_ct, rec]))
            pt_rows.append((pt, rgrp, risk_val, f"{s_ct} ({lig}) -> {t_ct} ({rec})", sc_v))
pt_comm = pd.DataFrame(pt_rows, columns=["patient", "sc_RiskGroup", "GMSig4_Risk", "interaction", "comm_score"])
pt_comm.to_csv(os.path.join(RES_DIR, "scRNA_TNBC_PatientLevel_CellComm.csv"), index=False)

# Print Patient-level correlations & Top Low vs High Differences
pt_cors = []
for inter, idf in pt_comm.groupby("interaction"):
    rho, sp_p = stats.spearmanr(idf["GMSig4_Risk"], idf["comm_score"])
    lo_m = idf.loc[idf["sc_RiskGroup"]=="Low_GMSig4", "comm_score"].mean()
    hi_m = idf.loc[idf["sc_RiskGroup"]=="High_GMSig4", "comm_score"].mean()
    _, wp = stats.mannwhitneyu(idf.loc[idf["sc_RiskGroup"]=="Low_GMSig4", "comm_score"],
                               idf.loc[idf["sc_RiskGroup"]=="High_GMSig4", "comm_score"])
    pt_cors.append((inter, lo_m, hi_m, np.log2((hi_m+0.01)/(lo_m+0.01)), rho, sp_p, wp))
pt_cor_df = pd.DataFrame(pt_cors, columns=["interaction", "mean_Low", "mean_High", "log2FC", "spearman_rho", "spearman_p", "wilcox_p"]).sort_values("spearman_rho")
pt_cor_df.to_csv(os.path.join(RES_DIR, "scRNA_TNBC_PatientLevel_CellComm_Summary.csv"), index=False)
print("\n=== Patient-Level (n=10) Cell-Cell Communication vs GMSig-4 Score ===")
print(pt_cor_df.to_string(index=False))

# Also run TME_niche sub-cluster Myeloid/T-cell/CAF focus interaction check (filtered to immune+tumor+CAF sub-niches)
sub_keep = ["CD8+ T", "CD4+ T", "Cycling T", "NK/NKT", "Macrophage", "DCs", "Monocyte", "iCAF", "myCAF", "Cancer Basal", "Cancer Cycling", "Endothelial"]
df_sub = df_ex[df_ex["TME_niche"].isin(sub_keep)].copy()
comm_sub_low  = compute_comm(df_sub[df_sub["sc_RiskGroup"]=="Low_GMSig4"],  "TME_niche", LR_PAIRS, n_perm=300, seed=303)
comm_sub_high = compute_comm(df_sub[df_sub["sc_RiskGroup"]=="High_GMSig4"], "TME_niche", LR_PAIRS, n_perm=300, seed=404)
comm_sub_m = comm_sub_low[["source", "target", "lr_pair", "pathway", "comm_score", "perm_p"]].rename(columns={"comm_score":"score_Low","perm_p":"p_Low"}).merge(
    comm_sub_high[["source", "target", "lr_pair", "pathway", "comm_score", "perm_p"]].rename(columns={"comm_score":"score_High","perm_p":"p_High"}),
    on=["source", "target", "lr_pair", "pathway"], how="outer"
).fillna({"score_Low":0.0, "score_High":0.0, "p_Low":1.0, "p_High":1.0})
comm_sub_m["log2FC_HighVsLow"] = np.log2((comm_sub_m["score_High"]+0.01)/(comm_sub_m["score_Low"]+0.01))
comm_sub_m.to_csv(os.path.join(RES_DIR, "scRNA_TNBC_SubNiche_CellComm_HighVsLow.csv"), index=False)

# ============================================================
# Render Supplementary Figure S3: scRNA Cell-Cell Communication
# ============================================================
fig = plt.figure(figsize=(18.4, 13.0), dpi=300)
gs = gridspec.GridSpec(2, 2, figure=fig, height_ratios=[1.10, 1.02], width_ratios=[1.0, 1.25], hspace=0.38, wspace=0.66)

# Panel A: Major Cell-Type Interaction Strength Delta (All permutation-active LR pairs p<0.05 in either group)
ax_a = fig.add_subplot(gs[0, 0])
active_comp = comp[(comp["p_Low"] < 0.05) | (comp["p_High"] < 0.05) | (comp["perm_fdr"] < 0.05)]
net_piv = active_comp.groupby(["source", "target"])["diff_HighMinusLow"].mean().unstack(fill_value=0)
net_cts = [c for c in ["T-cells", "Myeloid", "B-cells", "Plasmablasts", "CAFs", "Endothelial", "PVL", "Cancer Epithelial", "Normal Epithelial"] if c in net_piv.index or c in net_piv.columns]
net_piv = net_piv.reindex(index=net_cts, columns=net_cts, fill_value=0.0)
im_a = ax_a.imshow(net_piv.values, cmap="RdBu_r", vmin=-0.22, vmax=0.22)
ax_a.set_xticks(range(len(net_cts))); ax_a.set_xticklabels(net_cts, rotation=38, ha="right", fontsize=7.6)
ax_a.set_yticks(range(len(net_cts))); ax_a.set_yticklabels(net_cts, fontsize=7.8, fontweight="bold")
for r in range(len(net_cts)):
    for c in range(len(net_cts)):
        v = net_piv.values[r, c]
        if abs(v) >= 0.015:
            ax_a.text(c, r, f"{v:+.2f}", ha="center", va="center", fontsize=6.4, fontweight="bold",
                      color="white" if abs(v) > 0.13 else "#111")
cb_a = plt.colorbar(im_a, ax=ax_a, orientation="horizontal", fraction=0.042, pad=0.22)
cb_a.set_label("Mean Comm Score Delta (High - Low)", fontsize=7.5)
ax_a.set_xlabel("Target (Receptor Cell Type)")
ax_a.set_ylabel("Source (Ligand Cell Type)")
ax_a.set_title("A. Differential Cell-Cell Communication Matrix (High vs Low GMSig-4)", loc="left", fontweight="bold", fontsize=9.4)

# Panel B: Sub-Niche CXCL9/10-CXCR3, IFNG-IFNGR, CD274-PDCD1, SPP1/TGFB1 Bubble Plot (Low vs High GMSig-4)
ax_b = fig.add_subplot(gs[0, 1])
focus_pairs = [
    ("Macrophage", "CD8+ T", "CXCL10->CXCR3"),
    ("DCs",        "CD8+ T", "CXCL10->CXCR3"),
    ("Macrophage", "CD8+ T", "CXCL9->CXCR3"),
    ("Macrophage", "CD4+ T", "CXCL10->CXCR3"),
    ("CD8+ T",     "Macrophage", "IFNG->IFNGR1"),
    ("CD8+ T",     "DCs",        "IFNG->IFNGR1"),
    ("Macrophage", "CD8+ T", "CD274->PDCD1"),
    ("DCs",        "CD8+ T", "LGALS9->HAVCR2"),
    ("Macrophage", "CD8+ T", "CXCL16->CXCR6"),
    ("CD8+ T",     "Macrophage", "CCL5->CCR5"),
    ("Macrophage", "Cancer Basal", "SPP1->CD44"),
    ("Cancer Basal","Macrophage", "MIF->CD74"),
    ("iCAF",       "CD8+ T", "CXCL12->CXCR4"),
    ("myCAF",      "CD8+ T", "COL1A1->CD44"),
]
sub_rows = []
for s, t, lr in focus_pairs:
    m = comm_sub_m[(comm_sub_m["source"]==s) & (comm_sub_m["target"]==t) & (comm_sub_m["lr_pair"]==lr)]
    if len(m) > 0:
        row = m.iloc[0]
        lbl = f"{s} -> {t} | {lr}"
        sub_rows.append((lbl, "Low GMSig-4\n(Immune-Hot)", row["score_Low"], row["p_Low"]))
        sub_rows.append((lbl, "High GMSig-4\n(Immune-Cold)", row["score_High"], row["p_High"]))
b_df = pd.DataFrame(sub_rows, columns=["pair_lbl", "group", "score", "pval"])
pair_ord = list(dict.fromkeys([f"{s} -> {t} | {lr}" for s, t, lr in focus_pairs]))[::-1]
b_df["y_i"] = b_df["pair_lbl"].map({p: i for i, p in enumerate(pair_ord)})
b_df["x_i"] = b_df["group"].map({"Low GMSig-4\n(Immune-Hot)": 0, "High GMSig-4\n(Immune-Cold)": 1})

sizes = np.clip(-np.log10(np.clip(b_df["pval"].values, 1e-4, 1.0)) * 55 + 25, 25, 220)
sc_b = ax_b.scatter(b_df["x_i"], b_df["y_i"], s=sizes, c=b_df["score"], cmap="YlOrRd", vmin=0, vmax=max(b_df["score"].max(), 1.2), edgecolors="#222", lw=0.6)
for _, rr in b_df.iterrows():
    ax_b.text(rr["x_i"] + (0.14 if rr["x_i"]==1 else -0.14), rr["y_i"], f"{rr['score']:.2f}",
              ha="left" if rr["x_i"]==1 else "right", va="center", fontsize=7.0, fontweight="bold")
ax_b.set_xlim(-0.55, 1.55)
ax_b.set_xticks([0, 1]); ax_b.set_xticklabels(["Low GMSig-4 (n=5)", "High GMSig-4 (n=5)"], fontsize=8.2, fontweight="bold")
ax_b.set_yticks(range(len(pair_ord))); ax_b.set_yticklabels(pair_ord, fontsize=7.4)
ax_b.grid(color="#eeeeee", lw=0.5)
cb_b = plt.colorbar(sc_b, ax=ax_b, fraction=0.045, pad=0.02)
cb_b.set_label("Interaction Strength sqrt(L*R)", fontsize=7.5)
ax_b.set_title("B. Myeloid-T Cell & CAF Niche Communication (Low vs High GMSig-4)", loc="left", fontweight="bold", fontsize=9.4)

# Panel C: Patient-Level Correlation (n=10 TNBC Patients) with GMSig-4 Score
ax_c = fig.add_subplot(gs[1, 0])
pt_cor_df_s = pt_cor_df.sort_values("spearman_rho")
cols_c = ["#3C5488" if r < 0 else "#E64B35" for r in pt_cor_df_s["spearman_rho"]]
y_c = np.arange(len(pt_cor_df_s))
ax_c.barh(y_c, pt_cor_df_s["spearman_rho"], color=cols_c, height=0.66, edgecolor="white")
ax_c.axvline(0, color="#333333", ls="--", lw=0.8)
for idx_r, (_, rr) in enumerate(pt_cor_df_s.iterrows()):
    star = "***" if rr["spearman_p"]<0.001 else ("**" if rr["spearman_p"]<0.01 else ("*" if rr["spearman_p"]<0.05 else ""))
    ha_c = "right" if rr["spearman_rho"] < 0 else "left"
    off_c = -0.025 if rr["spearman_rho"] < 0 else 0.025
    ax_c.text(rr["spearman_rho"] + off_c, idx_r, f"rho={rr['spearman_rho']:+.2f}{star} (p={rr['spearman_p']:.3f})",
              va="center", ha=ha_c, fontsize=7.0, fontweight="bold" if star else "normal")
ax_c.set_yticks(y_c); ax_c.set_yticklabels(pt_cor_df_s["interaction"], fontsize=7.4)
ax_c.set_xlim(-1.42, 1.32)
ax_c.set_xlabel("Spearman rho vs Single-Cell Pseudobulk GMSig-4 Risk (n=10 TNBC)")
ax_c.grid(axis="x", color="#eeeeee", lw=0.5)
ax_c.set_title("C. Patient-Level Communication Rewiring Along GMSig-4 Axis", loc="left", fontweight="bold", fontsize=9.4)

# Panel D: Patient Scatter Plots of Key Immune Axes (Myeloid CXCL10->CXCR3 & T-cell IFNG->IFNGR1)
inner_d = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=gs[1, 1], wspace=0.28)
scatter_axes = [
    ("Myeloid (CXCL10) -> T-cells (CXCR3)", "Myeloid CXCL10 -> T CXCR3"),
    ("T-cells (IFNG) -> Myeloid (IFNGR1)",  "T-cell IFNG -> Myeloid IFNGR1"),
]
for k, (inter_key, stitle) in enumerate(scatter_axes):
    ax_s = fig.add_subplot(inner_d[k])
    sub_p = pt_comm[pt_comm["interaction"] == inter_key]
    for grp, col in [("Low_GMSig4", "#00A087"), ("High_GMSig4", "#E64B35")]:
        g_df = sub_p[sub_p["sc_RiskGroup"] == grp]
        ax_s.scatter(g_df["GMSig4_Risk"], g_df["comm_score"], s=52, c=col, edgecolor="black", lw=0.7,
                     label=grp.replace("_", " "))
    for _, rr in sub_p.iterrows():
        dy = -0.028 if rr["patient"] == "CID3946" else 0.015
        dx = -0.35 if rr["patient"] in ["CID3946", "CID44041"] else 0.07
        ax_s.text(rr["GMSig4_Risk"] + dx, rr["comm_score"] + dy, rr["patient"], fontsize=6.3, color="#333")
    # Linear fit
    z = np.polyfit(sub_p["GMSig4_Risk"], sub_p["comm_score"], 1)
    xs_line = np.linspace(sub_p["GMSig4_Risk"].min() - 0.15, sub_p["GMSig4_Risk"].max() + 0.15, 50)
    ax_s.plot(xs_line, np.polyval(z, xs_line), ls="--", color="#3C5488", lw=1.2)
    rho, sp_p = stats.spearmanr(sub_p["GMSig4_Risk"], sub_p["comm_score"])
    ax_s.text(0.96, 0.94, f"rho = {rho:.2f}\np = {sp_p:.4f}", transform=ax_s.transAxes,
              ha="right", va="top", fontsize=7.8, fontweight="bold",
              bbox=dict(boxstyle="round,pad=0.22", fc="white", ec="#cccccc", alpha=0.9))
    ax_s.set_xlim(sub_p["GMSig4_Risk"].min() - 0.4, sub_p["GMSig4_Risk"].max() + 0.75)
    ax_s.set_xlabel("Patient GMSig-4 Risk Score", fontsize=8.0)
    ax_s.set_ylabel("Interaction Strength sqrt(L*R)", fontsize=8.0)
    ax_s.set_title(stitle, fontsize=8.8, fontweight="bold")
    ax_s.grid(color="#eeeeee", lw=0.5)
    if k == 0:
        ax_s.legend(loc="lower left", fontsize=7.2, frameon=True)
        ax_s.text(0.0, 1.12, "D. Patient-Level Loss of Myeloid-T Cell Axis with Rising GMSig-4",
                  transform=ax_s.transAxes, fontsize=9.4, fontweight="bold", ha="left")

fig.suptitle("Supplementary Figure S3. Single-Cell Cell-Cell Communication Rewiring in TNBC (GSE176078, 42,512 Cells, 10 Patients)",
             fontsize=12.0, fontweight="bold", x=0.04, ha="left", y=0.975)

fig.savefig(os.path.join(SUP_DIR, "FigureS3_scRNA_CellComm_Mechanism.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(SUP_DIR, "FigureS3_scRNA_CellComm_Mechanism.png"), dpi=300, bbox_inches="tight")
print("SAVED FigureS3_scRNA_CellComm_Mechanism.{pdf,png} successfully!")
