"""Publication-Grade Figure 7: Orthogonal Tissue Validation — Patient-Level Spatial
Transcriptomics (GSE210616, N=22 TNBC, 43 Visium sections) & Clinical Mass-Spec
Proteomics (CPTAC-BRCA, n=28 TNBC, 0 overlap with TCGA)."""

import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
from scipy.stats import spearmanr

for fpath in [
    "/mnt/c/Windows/Fonts/arial.ttf",
    "/mnt/c/Windows/Fonts/arialbd.ttf",
    "/mnt/c/Windows/Fonts/ariali.ttf",
]:
    if os.path.exists(fpath):
        fm.fontManager.addfont(fpath)

plt.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.edgecolor": "#1E293B",
    "axes.linewidth": 0.65,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.labelsize": 8.6,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

ROOT   = "/mnt/d/Bioinformatics/TNBC_Microbiome"
SP_DIR = os.path.join(ROOT, "02_results/JTM_upgrade/spatial_gse210616")
CP_DIR = os.path.join(ROOT, "02_results/JTM_upgrade/cptac")
FIG_DIR = os.path.join(ROOT, "03_figures")

INK   = "#1E293B"
MUTED = "#64748B"
LOW   = "#1B6B93"
HIGH  = "#C0392B"
NAVY  = "#1B365D"
EMER  = "#1F8A70"
AMBER = "#D9822B"

pt_sp   = pd.read_csv(os.path.join(SP_DIR, "patient_level_spatial_metrics.csv"))
st_sp   = pd.read_csv(os.path.join(SP_DIR, "spatial_statistical_summary.csv"))
rep_df  = pd.read_csv(os.path.join(SP_DIR, "representative_slides.csv"))
cp_pt   = pd.read_csv(os.path.join(CP_DIR, "cptac_tnbc_patient_table.csv"))
cp_sum  = pd.read_csv(os.path.join(CP_DIR, "cptac_protein_concordance_summary.csv"))

low_r  = rep_df[rep_df["role"] == "Low_GMSig4_Rep"].iloc[0]
high_r = rep_df[rep_df["role"] == "High_GMSig4_Rep"].iloc[0]
sp_low  = pd.read_csv(os.path.join(SP_DIR, f"spots_Low_GMSig4_Rep_{low_r['patient_id']}_{low_r['gsm']}.csv"))
sp_high = pd.read_csv(os.path.join(SP_DIR, f"spots_High_GMSig4_Rep_{high_r['patient_id']}_{high_r['gsm']}.csv"))

fig = plt.figure(figsize=(15.8, 11.2), dpi=300, facecolor="white")
outer = gridspec.GridSpec(2, 1, figure=fig, height_ratios=[1.05, 1.0], hspace=0.40, top=0.92, bottom=0.06, left=0.05, right=0.98)
row1  = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[1.05, 1.25], wspace=0.22)
row2  = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[1], width_ratios=[1.05, 1.15], wspace=0.46)

def add_hdr(ax, letter, title, sub="", y=1.18):
    ax.text(-0.02, y, f"{letter}   {title}", transform=ax.transAxes,
            fontsize=10.0, fontweight="bold", color=NAVY, va="bottom", ha="left")
    if sub:
        ax.text(-0.02, y - 0.065, sub, transform=ax.transAxes,
                fontsize=7.5, color=MUTED, va="bottom", ha="left")

# ==============================================================================
# PANEL A: Representative Visium Tissue Maps (Low vs High GMSig-4, CD8 Lineage)
# ==============================================================================
gs_a = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=row1[0], wspace=0.08, hspace=0.28)
comp_cols = {"Tumor_Nest": "#E2A3A0", "Peritumoral_Margin": "#F4D06F", "Distal_Stroma": "#CBD5E1"}

for ci, (sdf, rinfo, title_tag, tcol) in enumerate([
    (sp_low,  low_r,  f"Low GMSig-4 ({low_r['patient_id']}, score={low_r['GMSig4']:.2f})", LOW),
    (sp_high, high_r, f"High GMSig-4 ({high_r['patient_id']}, score={high_r['GMSig4']:+.2f})", HIGH),
]):
    ax_top = fig.add_subplot(gs_a[0, ci])
    for cname, cval in [("Distal_Stroma", "#CBD5E1"), ("Peritumoral_Margin", "#F4C467"), ("Tumor_Nest", "#D97A72")]:
        sub = sdf[sdf["compartment"] == cname]
        ax_top.scatter(sub["array_col"], -sub["array_row"], s=5.5, c=cval, lw=0, label=cname.replace("_", " "))
    ax_top.set_title(title_tag, fontsize=8.3, fontweight="bold", color=tcol, pad=4)
    ax_top.set_aspect("equal"); ax_top.axis("off")
    if ci == 0:
        ax_top.legend(loc="lower left", fontsize=6.4, frameon=True, facecolor="white", framealpha=0.9,
                      edgecolor="#CBD5E0", markerscale=1.7, handletextpad=0.3, borderpad=0.3)
        fig.text(0.045, 0.962, "A   GSE210616 Visium Architecture & GMSig-4-Free CD8+ T-Cell Maps",
                 fontsize=10.2, fontweight="bold", color=NAVY, ha="left", va="bottom")
        fig.text(0.045, 0.943, "Compartments defined strictly without IFNG/IDO1/CXCL10/ENO1 (Epi vs Str; <=2 hex steps)",
                 fontsize=7.6, color=MUTED, ha="left", va="bottom")

    ax_bot = fig.add_subplot(gs_a[1, ci])
    vmax_cd8 = np.percentile(sp_low["cd8_score"].values, 98)
    sc = ax_bot.scatter(sdf["array_col"], -sdf["array_row"], s=5.8, c=sdf["cd8_score"],
                        cmap="viridis", vmin=0, vmax=max(vmax_cd8, 0.5), lw=0)
    ax_bot.set_title(f"CD8 T-cell (CD3D/E, CD8A/B): overall={rinfo['C_overall']:.3f}", fontsize=7.6, fontweight="bold", color=INK, pad=3)
    ax_bot.set_aspect("equal"); ax_bot.axis("off")
    if ci == 1:
        cax = ax_bot.inset_axes([1.02, 0.12, 0.045, 0.76])
        cb = fig.colorbar(sc, cax=cax)
        cb.set_label("CD8 log1p(CP10K)", fontsize=6.8, fontweight="bold", color=INK)
        cb.ax.tick_params(labelsize=6.5)

# ==============================================================================
# PANEL B: Patient-Level Spatial Quantitative Tests (N=22 Independent TNBC)
# ==============================================================================
gs_b = gridspec.GridSpecFromSubplotSpec(2, 2, subplot_spec=row1[1], wspace=0.28, hspace=0.40)
sp_panels = [
    ("C_overall", "Overall CD8+ T-Cell Score",        LOW,   "Pan-tissue immune-desert"),
    ("C_core",    "Tumor-Nest CD8+ Score (C_core)",   EMER,  "Depleted inside tumor nest"),
    ("C_margin",  "Peritumoral Margin CD8+ (C_margin)", AMBER, "Depleted at stromal margin"),
    ("CD8_PI",    "Core/Margin Ratio (CD8-PI)",       MUTED, "No selective margin trapping"),
]
for idx_p, (mcol, ylabel, col_p, tag_p) in enumerate(sp_panels):
    ax = fig.add_subplot(gs_b[idx_p // 2, idx_p % 2])
    r_all = st_sp[(st_sp["subset"] == "All_TNBC_n22") & (st_sp["score"] == "GMSig4") & (st_sp["metric"] == mcol)].iloc[0]
    r_cn  = st_sp[(st_sp["subset"] == "Chemo_Naive_n13") & (st_sp["score"] == "GMSig4") & (st_sp["metric"] == mcol)].iloc[0]

    for nac_lbl, mk, fc in [("Chemo_Naive", "o", col_p), ("NAC", "^", "white")]:
        sub = pt_sp[pt_sp["nac_status"] == nac_lbl]
        ax.scatter(sub["GMSig4"], sub[mcol], s=36, marker=mk, facecolor=fc, edgecolor=col_p, lw=1.1,
                   alpha=0.90, label=f"{nac_lbl.replace('_','-')} (n={len(sub)})", zorder=4)
    # Linear trendline for visual guide
    valid = pt_sp[["GMSig4", mcol]].dropna()
    z_fit = np.polyfit(valid["GMSig4"], valid[mcol], 1)
    xs = np.linspace(valid["GMSig4"].min(), valid["GMSig4"].max(), 50)
    ax.plot(xs, np.polyval(z_fit, xs), ls="--", color=col_p, lw=1.2, alpha=0.8)
    if mcol == "CD8_PI":
        ax.axhline(0, ls=":", color="#94A3B8", lw=0.8)

    ax.text(0.97, 0.96, f"All (N=22): rho={r_all['spearman_rho']:+.3f}, P={r_all['spearman_p']:.4f}\n"
                        f"Naive (N=13): rho={r_cn['spearman_rho']:+.3f}, P={r_cn['spearman_p']:.3f}",
            transform=ax.transAxes, ha="right", va="top", fontsize=6.5, fontweight="bold", color=INK,
            bbox=dict(boxstyle="round,pad=0.20", fc="white", ec="#CBD5E0", alpha=0.92, lw=0.5))
    ax.set_xlabel("Patient Pseudobulk GMSig-4", fontsize=7.6, fontweight="bold")
    ax.set_ylabel(ylabel, fontsize=7.5, fontweight="bold")
    ax.set_title(tag_p, fontsize=7.9, fontweight="bold", color=NAVY, pad=3)
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)
    if idx_p == 0:
        ax.legend(loc="lower left", fontsize=6.3, frameon=True, facecolor="white", edgecolor="#CBD5E0")
        fig.text(0.50, 0.962, "B   Patient-Level Spatial Analysis (GSE210616: N = 22 Patients, 43 Sections)",
                 fontsize=10.2, fontweight="bold", color=NAVY, ha="left", va="bottom")
        fig.text(0.50, 0.943, "Supports pan-compartment immune-desert depletion (C_core & C_margin both down) rather than margin trapping",
                 fontsize=7.6, color=MUTED, ha="left", va="bottom")

# ==============================================================================
# PANEL C: CPTAC RNA-to-Protein Concordance (n=28 TNBC, 0% TCGA Overlap)
# ==============================================================================
gs_c = gridspec.GridSpecFromSubplotSpec(1, 3, subplot_spec=row2[0], wspace=0.30)
conc_genes = [("ENO1", AMBER), ("IDO1", LOW), ("CXCL10", EMER)]
for gi, (g, col_g) in enumerate(conc_genes):
    ax = fig.add_subplot(gs_c[0, gi])
    sub = cp_pt[[f"RNA_{g}", f"PROT_{g}", "GMSig4_Group"]].dropna()
    rc = cp_sum[(cp_sum["cohort"] == "TNBC_All_n28") & (cp_sum["analysis_type"] == "RNA_protein_concordance") & (cp_sum["feature"] == g)].iloc[0]
    for grp_lbl, gc in [("Low", LOW), ("High", HIGH)]:
        sg = sub[sub["GMSig4_Group"] == grp_lbl]
        ax.scatter(sg[f"RNA_{g}"], sg[f"PROT_{g}"], s=30, c=gc, edgecolor="white", lw=0.5, alpha=0.9, label=f"{grp_lbl} GMSig-4")
    zf = np.polyfit(sub[f"RNA_{g}"], sub[f"PROT_{g}"], 1)
    xs = np.linspace(sub[f"RNA_{g}"].min(), sub[f"RNA_{g}"].max(), 40)
    ax.plot(xs, np.polyval(zf, xs), ls="--", color=NAVY, lw=1.1)
    p_str = f"{rc['spearman_p']:.1e}" if rc["spearman_p"] < 0.001 else f"{rc['spearman_p']:.4f}"
    ax.text(0.05, 0.95, f"n={int(rc['n'])}/28\nrho={rc['spearman_rho']:+.3f}\nP={p_str}",
            transform=ax.transAxes, ha="left", va="top", fontsize=6.9, fontweight="bold", color=INK,
            bbox=dict(boxstyle="round,pad=0.2", fc="white", ec="#CBD5E0", alpha=0.9, lw=0.5))
    ax.set_title(f"{g} (RNA vs MS)", fontsize=8.2, fontweight="bold", fontstyle="italic", color=NAVY)
    ax.set_xlabel("RNA log2 FPKM", fontsize=7.5, fontweight="bold")
    if gi == 0:
        ax.set_ylabel("MS Protein Abundance", fontsize=7.6, fontweight="bold")
        fig.text(0.045, 0.468, "C   CPTAC TNBC RNA-to-Protein Concordance (n = 28; 0% TCGA Overlap)",
                 fontsize=10.2, fontweight="bold", color=NAVY, ha="left", va="bottom")
        fig.text(0.045, 0.449, "Non-imputed MS proteomics; IFNG & CD274 proteins absent (0/28) and not imputed",
                 fontsize=7.6, color=MUTED, ha="left", va="bottom")
    if gi == 2:
        ax.legend(loc="lower right", fontsize=6.2, frameon=True, facecolor="white", edgecolor="#CBD5E0")
    ax.spines["top"].set_visible(False); ax.spines["right"].set_visible(False)

# ==============================================================================
# PANEL D: CPTAC Frozen RNA GMSig-4 vs Orthogonal Mass-Spec Proteins & TME
# ==============================================================================
ax_d = fig.add_subplot(row2[1])
feat_order = [
    ("Immune_4Prot_Avail_MeanZ",            "Immune 4-Prot (n=28)",        EMER),
    ("STAT1_protein",                       "STAT1 protein (n=28)",        EMER),
    ("CIBERSORT_ABSOLUTE_SCORE",            "CIBERSORT Immune (n=28)",     LOW),
    ("ESTIMATE_IMMUNE_SCORE",               "ESTIMATE Immune (n=28)",      LOW),
    ("CXCL10_protein",                      "CXCL10 protein (n=21)",       LOW),
    ("IDO1_protein",                        "IDO1 protein (n=23)",         LOW),
    ("GZMA_protein",                        "GZMA protein (n=28)",         EMER),
    ("CD8A_protein",                        "CD8A protein (n=27)",         EMER),
    ("PRF1_protein",                        "PRF1 protein (n=22)",         EMER),
    ("CD3_TILS_COUNTS",                     "Pathology CD3+ TILs (n=13)",  AMBER),
    ("ESTIMATE_STROMAL_SCORE",              "ESTIMATE Stromal (n=28)",     MUTED),
    ("ENO1_protein",                        "ENO1 protein (n=28)",         MUTED),
    ("Glycolysis_4Prot_Proxy_MeanZ_noENO1", "Glycolysis 4-Prot (n=28)",    MUTED),
]
sub_d = cp_sum[(cp_sum["cohort"] == "TNBC_All_n28") & (cp_sum["analysis_type"] == "GMSig4_association")].set_index("feature")
yd = np.arange(len(feat_order))
rhos = [float(sub_d.loc[f, "spearman_rho"]) for f, _, _ in feat_order]
pvs  = [float(sub_d.loc[f, "spearman_p"]) for f, _, _ in feat_order]
cols = [c for _, _, c in feat_order]

ax_d.axvline(0, color=INK, lw=0.7)
ax_d.barh(yd, rhos, height=0.66, color=cols, alpha=0.86, edgecolor="white", lw=0.6)
for yi, r, p in zip(yd, rhos, pvs):
    p_txt = f"P={p:.1e}" if p < 0.001 else f"P={p:.3f}"
    if r < 0:
        ax_d.text(r - 0.02, yi, f"{r:+.2f} ({p_txt})", va="center", ha="right", fontsize=6.7, fontweight="bold" if p < 0.05 else "normal", color=INK)
    else:
        ax_d.text(r + 0.02, yi, f"{r:+.2f} ({p_txt})", va="center", ha="left", fontsize=6.7, fontweight="bold" if p < 0.05 else "normal", color=MUTED)

ax_d.set_yticks(yd)
ax_d.set_yticklabels([lbl for _, lbl, _ in feat_order], fontsize=7.5, fontweight="bold", color=INK)
ax_d.invert_yaxis()
ax_d.set_xlim(-1.28, 0.42)
ax_d.set_xlabel("Spearman Correlation (rho) with Frozen RNA GMSig-4 Score", fontsize=8.0, fontweight="bold")
ax_d.spines["top"].set_visible(False); ax_d.spines["right"].set_visible(False)
fig.text(0.54, 0.468, "D   CPTAC Proteomic & TME Correlates of Frozen GMSig-4 (n = 28 TNBC)",
         fontsize=10.2, fontweight="bold", color=NAVY, ha="left", va="bottom")
fig.text(0.54, 0.449, "Strong depletion of cytotoxic/IFN immune proteins (P < 1e-8); bulk stromal & glycolysis proteins null",
         fontsize=7.6, color=MUTED, ha="left", va="bottom")

fig.savefig(os.path.join(FIG_DIR, "Figure7_Spatial_Proteomics.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG_DIR, "Figure7_Spatial_Proteomics.pdf"), bbox_inches="tight", facecolor="white")
print("SAVED Figure7_Spatial_Proteomics.png and .pdf successfully!")
