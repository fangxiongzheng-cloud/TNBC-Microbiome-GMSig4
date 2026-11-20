import os
import pandas as pd
import numpy as np
import matplotlib

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
import matplotlib.patches as patches
import matplotlib.patheffects as pe
from matplotlib.colors import LinearSegmentedColormap
from sklearn.metrics import roc_curve, auc

# Use locally installed Arial faces when an explicit font directory is supplied.
arial_dir = os.environ.get("ARIAL_FONT_DIR")
if arial_dir:
    for filename in ("arial.ttf", "arialbd.ttf", "ariali.ttf", "arialbi.ttf"):
        font_path = os.path.join(arial_dir, filename)
        if os.path.isfile(font_path):
            fm.fontManager.addfont(font_path)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.edgecolor": "#2D3748",
    "axes.linewidth": 0.65,
    "xtick.labelsize": 7.8,
    "ytick.labelsize": 7.8,
    "axes.labelsize": 8.4,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

RES_DIR = os.path.join(project_root, "02_results")
FIG_DIR = os.path.join(project_root, "03_figures")
os.makedirs(FIG_DIR, exist_ok=True)

INK   = "#1A202C"
MUTED = "#4A5568"
TEAL  = "#2B7A9E"
EMER  = "#1F8A70"
NAVY  = "#1E3A5F"
CORAL = "#C84B31"
CRIM  = "#8B1E1E"

ici_sum   = pd.read_csv(os.path.join(RES_DIR, "GMSig4_Immunotherapy_pCR_Validation_Summary.csv"))
ici_exact = pd.read_csv(os.path.join(RES_DIR, "GMSig4_ICI_Subgroup_ExactScores.csv"))
durva     = pd.read_csv(os.path.join(RES_DIR, "GMSig4_PatientScores_GSE173839_Durva.csv"))
umap_df   = pd.read_csv(os.path.join(RES_DIR, "scRNA_TNBC_UMAP_metadata_keygenes.csv"))
sc_mean   = pd.read_csv(os.path.join(RES_DIR, "scRNA_TNBC_gutMGene_celltype_major_mean.csv"), index_col=0)
sc_pct    = pd.read_csv(os.path.join(RES_DIR, "scRNA_TNBC_gutMGene_celltype_major_pct.csv"), index_col=0)

def add_panel_header(ax, letter, title, y_axes=1.06, x_letter_pts=0):
    ax.annotate(letter, xy=(0, y_axes), xycoords=ax.transAxes, xytext=(x_letter_pts, 0),
                textcoords="offset points", fontsize=14, fontweight="bold", color="#111111", va="bottom", ha="left")
    ax.annotate(title, xy=(0, y_axes), xycoords=ax.transAxes, xytext=(x_letter_pts + 17, 1.2),
                textcoords="offset points", fontsize=9.8, fontweight="bold", color=INK, va="bottom", ha="left")

fig = plt.figure(figsize=(16.8, 12.2), dpi=300, facecolor="white")
outer = gridspec.GridSpec(2, 3, figure=fig, height_ratios=[0.95, 1.18], width_ratios=[1.05, 1.0, 1.08],
                          hspace=0.35, wspace=0.25)

ici_map = ici_sum.set_index("cohort")

# ==============================================================================
# PANEL A: Immunotherapy Specificity Raincloud (I-SPY2 GSE194040)
# ==============================================================================
ax5a = fig.add_subplot(outer[0, 0])
arm_defs = [
    ("Pembro\nTNBC\n(n=29)",          "Pembro_TNBC",          "GSE194040_Pembro_TNBC"),
    ("Pembro\nAll\n(n=69)",           "Pembro_All",           "GSE194040_Pembro_All"),
    ("Pembro\nHR+HER2\u2212\n(n=40)", "Pembro_HRp",           "GSE194040_Pembro_HRpHER2n"),
    ("Paclitaxel\nCtrl TNBC\n(n=85)", "Paclitaxel_Ctrl_TNBC", "GSE194040_Paclitaxel_Ctrl_TNBC"),
]
xticks, xlabels = [], []
rng = np.random.default_rng(42)

for idx, (lbl, sub_key, coh_key) in enumerate(arm_defs):
    sub = ici_exact[ici_exact["subgroup"] == sub_key].copy()
    pos0 = idx * 2.5
    pos1 = pos0 + 0.92
    v0 = sub.loc[sub["pCR"] == 0, "Ridge_Score"].dropna().values
    v1 = sub.loc[sub["pCR"] == 1, "Ridge_Score"].dropna().values
    xticks.append((pos0 + pos1) / 2.0)
    xlabels.append(lbl)

    for v, pos, col, side in [(v0, pos0, CORAL, -1), (v1, pos1, TEAL, 1)]:
        # Half violin
        vp = ax5a.violinplot([v], positions=[pos], widths=0.72, showmeans=False, showmedians=False, showextrema=False)
        for pc in vp["bodies"]:
            verts = pc.get_paths()[0].vertices
            if side == -1:
                verts[:, 0] = np.clip(verts[:, 0], -np.inf, pos)
            else:
                verts[:, 0] = np.clip(verts[:, 0], pos, np.inf)
            pc.set_facecolor(col); pc.set_alpha(0.30); pc.set_edgecolor(col); pc.set_linewidth(0.65)
        ax5a.boxplot([v], positions=[pos], widths=0.22, patch_artist=True, showfliers=False,
                     boxprops=dict(facecolor=col, alpha=0.85, edgecolor=INK, lw=0.7),
                     medianprops=dict(color="white", lw=1.2),
                     whiskerprops=dict(color=INK, lw=0.7), capprops=dict(color=INK, lw=0.7))
        ax5a.scatter(pos + rng.uniform(-0.10, 0.10, len(v)), v, s=11, color=col, edgecolor="white", lw=0.35, alpha=0.85, zorder=4)

    wp = float(ici_map.loc[coh_key, "GMSig4_wilcox_p"])
    ymax = max(np.max(v0), np.max(v1))
    p_str = f"P={wp:.3f}" if wp >= 0.005 else f"P={wp:.2e}"
    sig_col = CRIM if wp < 0.05 else MUTED
    ax5a.plot([pos0, pos0, pos1, pos1], [ymax + 0.16, ymax + 0.27, ymax + 0.27, ymax + 0.16], color=INK, lw=0.75)
    ax5a.text((pos0 + pos1)/2.0, ymax + 0.32, p_str, ha="center", va="bottom",
              fontsize=7.1, fontweight="bold" if wp < 0.05 else "normal", color=sig_col)

ax5a.axhline(0, ls="--", color="#94A3B8", lw=0.7)
ax5a.set_xticks(xticks)
ax5a.set_xticklabels(xlabels, fontsize=7.4, fontweight="bold", color=INK)
ax5a.set_ylabel("GMSig-4 Risk Score", fontweight="bold", color=INK)
ax5a.set_ylim(-2.4, 2.75)
ax5a.spines["top"].set_visible(False); ax5a.spines["right"].set_visible(False)
ax5a.legend(handles=[patches.Patch(fc=TEAL, alpha=0.85, label="pCR (Responder)"),
                     patches.Patch(fc=CORAL, alpha=0.85, label="RD (Non-Responder)")],
            loc="lower left", fontsize=7.2, frameon=False, ncol=2)
add_panel_header(ax5a, "A", "Immunotherapy Specificity (I-SPY2 GSE194040)", y_axes=1.06, x_letter_pts=-6)

# ==============================================================================
# PANEL B: pCR Stratification w/ Fisher Exact P Brackets
# ==============================================================================
ax5b = fig.add_subplot(outer[0, 1])
sel_rows = [
    ("GSE194040_Pembro_TNBC",          "Pembro\nTNBC (n=29)"),
    ("GSE194040_Pembro_All",           "Pembro\nAll (n=69)"),
    ("GSE194040_Pembro_HRpHER2n",      "Pembro\nHR+ (n=40)"),
    ("GSE173839_DurvaOlap_TNBC",       "Durva+Olap\nTNBC (n=21)"),
    ("GSE194040_Paclitaxel_Ctrl_TNBC", "Paclitaxel\nCtrl (n=85)"),
]
sub_sum = ici_sum.set_index("cohort").loc[[c for c, _ in sel_rows]].copy()
x = np.arange(len(sel_rows))
w = 0.35
b_low  = ax5b.bar(x - w/2, sub_sum["pCR_rate_LowRisk"] * 100, w, color=TEAL, edgecolor="white", label="Low GMSig-4 Risk")
b_high = ax5b.bar(x + w/2, sub_sum["pCR_rate_HighRisk"] * 100, w, color=CORAL, edgecolor="white", label="High GMSig-4 Risk")

for i_r, (coh_k, _) in enumerate(sel_rows):
    yl = sub_sum.loc[coh_k, "pCR_rate_LowRisk"] * 100
    yh = sub_sum.loc[coh_k, "pCR_rate_HighRisk"] * 100
    fp = float(sub_sum.loc[coh_k, "Fisher_p"])
    ax5b.text(i_r - w/2, yl + 1.5, f"{yl:.1f}%", ha="center", va="bottom", fontsize=6.7, fontweight="bold", color=INK)
    ax5b.text(i_r + w/2, yh + 1.5, f"{yh:.1f}%", ha="center", va="bottom", fontsize=6.7, fontweight="bold", color=INK)
    top_y = max(yl, yh) + 7.5
    p_lbl = f"P={fp:.3f}" if fp < 0.05 else f"P={fp:.2f}"
    ax5b.plot([i_r - w/2, i_r - w/2, i_r + w/2, i_r + w/2], [top_y - 1.5, top_y, top_y, top_y - 1.5], color=INK, lw=0.7)
    ax5b.text(i_r, top_y + 1.0, p_lbl, ha="center", va="bottom", fontsize=6.7,
              fontweight="bold" if fp < 0.05 else "normal", color=CRIM if fp < 0.05 else MUTED)

ax5b.set_xticks(x)
ax5b.set_xticklabels([lbl for _, lbl in sel_rows], fontsize=7.2, fontweight="bold", color=INK)
ax5b.set_ylabel("pCR Rate (%)", fontweight="bold", color=INK)
ax5b.set_ylim(0, 102)
ax5b.spines["top"].set_visible(False); ax5b.spines["right"].set_visible(False)
ax5b.legend(loc="upper right", fontsize=7.2, frameon=False)
add_panel_header(ax5b, "B", "pCR Stratification Across ICI vs Control", y_axes=1.06, x_letter_pts=-6)

# ==============================================================================
# PANEL C: Immunotherapy pCR Prediction ROC
# ==============================================================================
ax5c = fig.add_subplot(outer[0, 2])
p_tn_df  = ici_exact[ici_exact["subgroup"] == "Pembro_TNBC"]
p_all_df = ici_exact[ici_exact["subgroup"] == "Pembro_All"]
d_tn_ex  = ici_exact[ici_exact["subgroup"] == "DurvaOlap_TNBC"]
c_tn_df  = ici_exact[ici_exact["subgroup"] == "Paclitaxel_Ctrl_TNBC"]
d_tn_df  = durva[durva["is_TNBC"] == True]

roc_specs = [
    ("GMSig-4 Pembro TNBC",    p_tn_df["pCR"],  -p_tn_df["Ridge_Score"],   CORAL,    "-",  2.3),
    ("UniBeta Pembro TNBC",    p_tn_df["pCR"],  -p_tn_df["UniBeta_Score"], CRIM,     "--", 1.7),
    ("GMSig-4 Pembro All",     p_all_df["pCR"], -p_all_df["Ridge_Score"],  EMER,     "-",  1.9),
    ("GMSig-4 Durva TNBC",     d_tn_ex["pCR"],  -d_tn_ex["Ridge_Score"],   NAVY,     "-",  1.8),
    ("CD274 (PD-L1) Durva",    d_tn_df["pCR"],   d_tn_df["PDL1"],          "#7F8C8D",":",  1.5),
    ("Paclitaxel Ctrl TNBC",   c_tn_df["pCR"],  -c_tn_df["Ridge_Score"],   "#B0BAC5","-.", 1.5),
]
ax5c.plot([0, 1], [0, 1], ls="--", color="#CBD5E0", lw=0.9)
for lbl, y_t, y_s, col, ls, lw in roc_specs:
    fpr, tpr, _ = roc_curve(y_t, y_s)
    a_val = auc(fpr, tpr)
    ax5c.plot(fpr, tpr, color=col, ls=ls, lw=lw, label=f"{lbl} ({a_val:.3f})")

ax5c.set_xlabel("False Positive Rate (1 \u2212 Specificity)", fontweight="bold", color=INK)
ax5c.set_ylabel("True Positive Rate (Sensitivity)", fontweight="bold", color=INK)
ax5c.set_xlim(-0.02, 1.02); ax5c.set_ylim(-0.02, 1.02)
ax5c.spines["top"].set_visible(False); ax5c.spines["right"].set_visible(False)
ax5c.legend(title="Cohort (AUC)", title_fontsize=7.2, loc="lower right", fontsize=6.9,
            frameon=True, facecolor="white", framealpha=0.92, edgecolor="#CBD5E0")
add_panel_header(ax5c, "C", "Immunotherapy pCR Prediction ROC", y_axes=1.06, x_letter_pts=-6)

# ==============================================================================
# PANEL D: Borderless 42,512-Cell TNBC scRNA-seq Atlas w/ Corner UMAP Arrows
# ==============================================================================
ax5d = fig.add_subplot(outer[1, 0])
ct_palette = {
    "T-cells":           "#C84B31",
    "Cancer Epithelial": "#1E3A5F",
    "Myeloid":           "#1F8A70",
    "Plasmablasts":      "#E8A87C",
    "CAFs":              "#7E8B9B",
    "B-cells":           "#2B7A9E",
    "Endothelial":       "#8C6D46",
    "PVL":               "#82C09A",
    "Normal Epithelial": "#D9822B",
}
ct_counts = umap_df["celltype_major"].value_counts()
shuf = umap_df.sample(frac=1.0, random_state=7)
for ct, col in ct_palette.items():
    sub = shuf[shuf["celltype_major"] == ct]
    ax5d.scatter(sub["UMAP_1"], sub["UMAP_2"], s=0.9, c=col, alpha=0.68, rasterized=True, edgecolors="none")
    cx, cy = sub["UMAP_1"].median(), sub["UMAP_2"].median()
    ax5d.text(cx, cy, f"{ct}\n({ct_counts[ct]:,})", fontsize=6.8, fontweight="bold", ha="center", va="center",
              color=INK, path_effects=[pe.withStroke(linewidth=2.4, foreground="white", alpha=0.92)])

ax5d.set_xticks([]); ax5d.set_yticks([])
for sp in ax5d.spines.values(): sp.set_visible(False)
# Sleek L-shaped UMAP corner arrow indicator
x0, y0 = -8.5, -12.0
ax5d.annotate("", xy=(x0 + 6.5, y0), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.1, mutation_scale=9))
ax5d.annotate("", xy=(x0, y0 + 5.5), xytext=(x0, y0), arrowprops=dict(arrowstyle="-|>", color=INK, lw=1.1, mutation_scale=9))
ax5d.text(x0 + 3.2, y0 - 1.3, "UMAP 1", ha="center", va="top", fontsize=6.8, fontweight="bold", color=INK)
ax5d.text(x0 - 1.2, y0 + 2.7, "UMAP 2", ha="right", va="center", rotation=90, fontsize=6.8, fontweight="bold", color=INK)
ax5d.set_xlim(-10.5, 27.5); ax5d.set_ylim(-14.5, 22.5)
add_panel_header(ax5d, "D", "TNBC scRNA-seq Atlas (42,512 Cells, n = 10)", y_axes=1.02, x_letter_pts=-6)

# ==============================================================================
# PANEL E: Single-Cell Niche Mapping (Fluorescence-Style Grey->Crimson + 1 Colorbar)
# ==============================================================================
gs_e = gridspec.GridSpecFromSubplotSpec(3, 2, subplot_spec=outer[1, 1], height_ratios=[1.0, 1.0, 0.14], hspace=0.20, wspace=0.08)
sig_genes = [
    ("IFNG",   "IFNG (T-cells 23.9%)"),
    ("IDO1",   "IDO1 (Myeloid 9.3%)"),
    ("CXCL10", "CXCL10 (Myeloid 33.7%)"),
    ("ENO1",   "ENO1 (Cancer Epi 90.5%)"),
]
expr_cmap = LinearSegmentedColormap.from_list("sc_expr", ["#EAECEF", "#F4A582", "#D6604D", "#99000D"])

last_sc = None
for k, (g, gtitle) in enumerate(sig_genes):
    ax_g = fig.add_subplot(gs_e[k // 2, k % 2])
    vals = umap_df[g].values
    vmax = np.percentile(vals, 99.2)
    vals_norm = np.clip(vals / max(vmax, 1e-6), 0, 1)
    order = np.argsort(vals_norm)
    last_sc = ax_g.scatter(umap_df["UMAP_1"].values[order], umap_df["UMAP_2"].values[order],
                           c=vals_norm[order], cmap=expr_cmap, s=0.65, vmin=0, vmax=1,
                           rasterized=True, edgecolors="none")
    ax_g.set_title(gtitle, fontsize=7.6, fontweight="bold", style="italic", color=INK, pad=2)
    ax_g.set_xticks([]); ax_g.set_yticks([])
    for sp in ax_g.spines.values(): sp.set_visible(False)
    if k == 0:
        add_panel_header(ax_g, "E", "Single-Cell Niche Mapping of GMSig-4", y_axes=1.16, x_letter_pts=-6)

ax_e_cb = fig.add_subplot(gs_e[2, :])
ax_e_cb.axis("off")
cax_e = ax_e_cb.inset_axes([0.18, 0.40, 0.64, 0.42])
cb_e = fig.colorbar(last_sc, cax=cax_e, orientation="horizontal")
cb_e.set_ticks([0, 0.5, 1.0])
cb_e.set_ticklabels(["Low (0)", "Mid", "High (99th pct)"])
cb_e.ax.tick_params(labelsize=6.8)
cb_e.outline.set_linewidth(0.45)

# ==============================================================================
# PANEL F: Gut-Host Causal & Signature Niche Dotplot w/ Lineage & Tier Strips
# ==============================================================================
gs_f = gridspec.GridSpecFromSubplotSpec(
    2, 2, subplot_spec=outer[1, 2],
    height_ratios=[0.045, 1.0], width_ratios=[0.045, 1.0],
    hspace=0.02, wspace=0.02
)
ax5f_top  = fig.add_subplot(gs_f[0, 1])
ax5f_left = fig.add_subplot(gs_f[1, 0])
ax5f      = fig.add_subplot(gs_f[1, 1])

dot_genes = ["ENO1", "LDHA", "MUC1", "CLDN7", "TCF7L2", "PEMT", "DAPK3", "MAGI3",
             "TGFBR1", "FFAR2", "CASP1", "IDO1", "CXCL10", "IFNG", "CXCR6", "FOXP3", "CCR6"]
ct_order  = ["Cancer Epithelial", "Normal Epithelial", "CAFs", "PVL", "Endothelial",
             "Myeloid", "T-cells", "B-cells", "Plasmablasts"]
ct_short  = ["Cancer Epi", "Normal Epi", "CAFs", "PVL", "Endothelial", "Myeloid", "T-cells", "B-cells", "Plasmablasts"]
ct_group  = ["Epithelial", "Epithelial", "Stromal", "Stromal", "Stromal", "Immune", "Immune", "Immune", "Immune"]
grp_cols  = {"Epithelial": NAVY, "Stromal": "#8C6D46", "Immune": EMER}

dot_genes = [g for g in dot_genes if g in sc_mean.index]
mean_sub = sc_mean.loc[dot_genes, ct_order]
pct_sub  = sc_pct.loc[dot_genes, ct_order] * 100.0
mean_z   = mean_sub.sub(mean_sub.mean(axis=1), axis=0).div(mean_sub.std(axis=1) + 1e-8, axis=0)

# Top Lineage Strip
for xi, cg in enumerate(ct_group):
    ax5f_top.add_patch(patches.Rectangle((xi - 0.5, 0), 1.0, 1.0, fc=grp_cols[cg], ec="white", lw=0.8))
ax5f_top.set_xlim(-0.5, len(ct_order) - 0.5); ax5f_top.set_ylim(0, 1); ax5f_top.axis("off")
add_panel_header(ax5f_top, "F", "Gut\u2013Host Causal & Signature Niche Dotplot", y_axes=1.55, x_letter_pts=-20)

# Left Gene Tier Strip (GMSig-4 vs MR/Docking Target)
gmsig_set = {"ENO1", "IDO1", "CXCL10", "IFNG"}
for yi, g in enumerate(reversed(dot_genes)):
    fc_t = CORAL if g in gmsig_set else "#94A3B8"
    ax5f_left.add_patch(patches.Rectangle((0, yi - 0.5), 1.0, 1.0, fc=fc_t, ec="white", lw=0.8))
ax5f_left.set_xlim(0, 1); ax5f_left.set_ylim(-0.5, len(dot_genes) - 0.5)
ax5f_left.set_yticks(range(len(dot_genes)))
ax5f_left.set_yticklabels(list(reversed(dot_genes)), fontsize=7.5, fontweight="bold", fontstyle="italic", color=INK)
ax5f_left.set_xticks([]); ax5f_left.tick_params(axis="y", length=0, pad=3)
for sp in ax5f_left.spines.values(): sp.set_visible(False)

xs, ys, sizes, cols = [], [], [], []
for yi, g in enumerate(reversed(dot_genes)):
    for xi, ct in enumerate(ct_order):
        xs.append(xi); ys.append(yi)
        sizes.append(max(pct_sub.loc[g, ct] * 1.45, 4.0))
        cols.append(mean_z.loc[g, ct])

dot_cmap = LinearSegmentedColormap.from_list("dot_div", ["#1D5B9B", "#F8FAFC", "#C0392B"])
sc_dot = ax5f.scatter(xs, ys, s=sizes, c=cols, cmap=dot_cmap, vmin=-1.6, vmax=1.8, edgecolors=INK, linewidths=0.35)
ax5f.set_xlim(-0.5, len(ct_order) - 0.5); ax5f.set_ylim(-0.5, len(dot_genes) - 0.5)
ax5f.set_xticks(range(len(ct_short)))
ax5f.set_xticklabels(ct_short, rotation=36, ha="right", fontsize=7.3, fontweight="bold", color=INK)
ax5f.set_yticks([])
ax5f.grid(color="#EDF2F7", lw=0.5)
for sp in ax5f.spines.values(): sp.set_color("#CBD5E0"); sp.set_linewidth(0.5)

cbar_f = plt.colorbar(sc_dot, ax=ax5f, fraction=0.042, pad=0.03)
cbar_f.set_label("Mean Expr (Row Z)", fontsize=7.3, fontweight="bold")
cbar_f.ax.tick_params(labelsize=6.8)
cbar_f.outline.set_linewidth(0.45)

for pct_ref in [10, 40, 80]:
    ax5f.scatter([], [], s=pct_ref * 1.45, c="#94A3B8", edgecolors=INK, lw=0.4, label=f"{pct_ref}%")
ax5f.legend(title="% Cells", title_fontsize=7.0, fontsize=6.8, loc="upper center",
            bbox_to_anchor=(0.48, -0.22), frameon=False, ncol=3, columnspacing=0.9, handletextpad=0.3)

fig.savefig(os.path.join(FIG_DIR, "Figure5.png"), dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(os.path.join(FIG_DIR, "Figure5.pdf"), bbox_inches="tight", facecolor="white")
print("SAVED Figure5.png and Figure5.pdf successfully!")
