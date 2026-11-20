import os
import io
import pandas as pd
import numpy as np
import matplotlib

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.font_manager as fm
import matplotlib.patches as patches
from matplotlib.colors import LinearSegmentedColormap
from scipy.cluster.hierarchy import linkage, dendrogram
from scipy.spatial.distance import pdist
from PIL import Image

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.Draw import rdMolDraw2D

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
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.labelsize": 8.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

RES_DIR    = os.path.join(project_root, "02_results")
FIG_DIR    = os.path.join(project_root, "03_figures")
ASSET_DIR = os.path.join(project_root, "01_data", "docking", "figure_assets")
POCKET_DIR = os.path.join(ASSET_DIR, "pockets3d")
BG_6D_PATH = os.path.join(ASSET_DIR, "panel6d_background.jpg")
if not os.path.isdir(POCKET_DIR) or not os.path.isfile(BG_6D_PATH):
    raise FileNotFoundError(
        "Figure 6 v2 requires local illustration assets. "
        "Run 22_fig6_docking_and_mechanism.py for the data-driven figure."
    )
os.makedirs(FIG_DIR, exist_ok=True)

INK   = "#1A202C"
MUTED = "#4A5568"
TEAL  = "#2B7A9E"
EMER  = "#1F8A70"
NAVY  = "#1E3A5F"
CORAL = "#C84B31"
OCHRE = "#C07A28"

dock_df = pd.read_csv(os.path.join(RES_DIR, "Docking_Vina_Affinities_Matrix.csv"))
prot_df = pd.read_csv(os.path.join(RES_DIR, "Docking_Protein_Manifest.csv"))

def render_rdkit_png(smiles, width=380, height=120):
    mol = Chem.MolFromSmiles(smiles)
    AllChem.Compute2DCoords(mol)
    drawer = rdMolDraw2D.MolDraw2DCairo(width, height)
    opts = drawer.drawOptions()
    opts.clearBackground = False
    opts.bondLineWidth = 2.0
    opts.padding = 0.06
    opts.minFontSize = 13
    drawer.DrawMolecule(mol)
    drawer.FinishDrawing()
    img = Image.open(io.BytesIO(drawer.GetDrawingText())).convert("RGBA")
    arr = np.array(img)
    white = (arr[:, :, 0] > 246) & (arr[:, :, 1] > 246) & (arr[:, :, 2] > 246)
    arr[white, 3] = 0
    return Image.fromarray(arr)

def add_panel_header(ax, letter, title, y_axes=1.05, x_letter_pts=0):
    ax.annotate(letter, xy=(0, y_axes), xycoords=ax.transAxes, xytext=(x_letter_pts, 0),
                textcoords="offset points", fontsize=14, fontweight="bold", color="#111111", va="bottom", ha="left")
    ax.annotate(title, xy=(0, y_axes), xycoords=ax.transAxes, xytext=(x_letter_pts + 17, 1.2),
                textcoords="offset points", fontsize=10.0, fontweight="bold", color=INK, va="bottom", ha="left")

# ==============================================================================
# Build Figure 6 Canvas
# ==============================================================================
fig = plt.figure(figsize=(15.6, 17.2), dpi=300, facecolor="white")
outer = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1.02, 0.92, 0.94], hspace=0.58)

top_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[1.38, 1.0], wspace=0.22)
mid_gs = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[1], wspace=0.14)

# ==============================================================================
# PANEL A: Clustered 11x11 Vina Matrix + Dendrograms + Strong-Pocket Margin
# ==============================================================================
gs_a = gridspec.GridSpecFromSubplotSpec(
    3, 4, subplot_spec=top_gs[0],
    height_ratios=[0.16, 0.06, 1.0],
    width_ratios=[0.14, 0.045, 1.0, 0.24],
    hspace=0.025, wspace=0.02
)
ax_a_cdend = fig.add_subplot(gs_a[0, 2])
ax_a_top   = fig.add_subplot(gs_a[1, 2])
ax_a_rdend = fig.add_subplot(gs_a[2, 0])
ax_a_left  = fig.add_subplot(gs_a[2, 1])
ax_a_mat   = fig.add_subplot(gs_a[2, 2])
ax_a_rbar  = fig.add_subplot(gs_a[2, 3])

mat_raw = dock_df.pivot(index="gene", columns="metabolite", values="affinity_kcal_mol")
ev_raw  = dock_df.fillna({"gutMGene_evidence": ""}).pivot(index="gene", columns="metabolite", values="gutMGene_evidence")

# Hierarchical clustering on rows & columns
row_lnk = linkage(pdist(mat_raw.values, metric="euclidean"), method="ward")
col_lnk = linkage(pdist(mat_raw.values.T, metric="euclidean"), method="ward")

d_row = dendrogram(row_lnk, orientation="left", ax=ax_a_rdend, no_labels=True,
                   color_threshold=0, above_threshold_color="#64748B")
d_col = dendrogram(col_lnk, orientation="top", ax=ax_a_cdend, no_labels=True,
                   color_threshold=0, above_threshold_color="#64748B")

for ax_d in (ax_a_rdend, ax_a_cdend):
    for col_line in ax_d.collections:
        col_line.set_linewidth(0.85)
    ax_d.axis("off")

# Note: imshow y=0 is top; matplotlib dendrogram leaves go bottom-to-top, so reverse leaves for top-to-bottom alignment
gene_order = [mat_raw.index[i] for i in d_row["leaves"][::-1]]
met_order  = [mat_raw.columns[i] for i in d_col["leaves"]]

mat_piv = mat_raw.loc[gene_order, met_order]
ev_piv  = ev_raw.loc[gene_order, met_order]

class_map = {
    "Lithocholic_acid": "Bile Acid",
    "Urolithin_A": "Polyphenol",
    "Caffeic_acid": "Polyphenol",
    "3-(4-Hydroxyphenyl)propionic_acid": "Polyphenol",
    "3-Indolepropionic_acid": "Indole",
    "Indole-3-carboxylic_acid": "Indole",
    "Indole": "Indole",
    "Nicotinic_acid": "B3 / Choline",
    "Butyrate": "SCFA",
    "Propionate": "SCFA",
    "Trimethylamine_N-oxide": "B3 / Choline",
}
class_colors = {
    "Bile Acid": "#8C6D46",
    "Polyphenol": "#1F8A70",
    "Indole": "#3C5488",
    "B3 / Choline": "#7E8B9B",
    "SCFA": "#2B7A9E",
}
met_disp = {
    "Lithocholic_acid": "Lithocholic acid",
    "Urolithin_A": "Urolithin A",
    "3-Indolepropionic_acid": "IPA",
    "Caffeic_acid": "Caffeic acid",
    "Indole-3-carboxylic_acid": "Indole-3-CA",
    "3-(4-Hydroxyphenyl)propionic_acid": "4-HPPA",
    "Indole": "Indole",
    "Nicotinic_acid": "Nicotinic acid",
    "Butyrate": "Butyrate",
    "Propionate": "Propionate",
    "Trimethylamine_N-oxide": "TMAO",
}

# Top metabolite scaffold color bar
for c_i, m in enumerate(met_order):
    ax_a_top.add_patch(patches.Rectangle((c_i - 0.5, 0), 1.0, 1.0, fc=class_colors[class_map[m]], ec="white", lw=1.0))
ax_a_top.set_xlim(-0.5, len(met_order) - 0.5); ax_a_top.set_ylim(0, 1)
ax_a_top.axis("off")

add_panel_header(ax_a_cdend, "A", "Clustered AutoDock Vina Affinity Matrix (11 Targets \u00D7 11 Metabolites)",
                 y_axes=1.12, x_letter_pts=-68)

# Left structure source strip & protein gene names aligned to left of rdend? Wait: put labels between rdend?
# Better: put gene labels cleanly on ax_a_left (between rdend and matrix) or right of rbar!
# Let's check: putting gene labels on ax_a_rbar y-axis (left) or ax_a_left: ax_a_rdend is at [2,0], ax_a_left is [2,1].
# If we put gene labels on ax_a_rdend left? No, dendrogram branches go left-to-right into ax_a_left!
# Putting italic gene labels on the RIGHT of ax_a_rbar? Or putting structure strip at [2,0], gene labels on ax_a_left?
# Wait: standard ComplexHeatmap puts: [Gene Labels] [Source Strip] [Heatmap] [Row Bar] [Right Dendrogram] OR [Left Dendrogram] [Source Strip] [Heatmap] [Right Gene Labels + Row Bar]!
# Putting gene names between Source Strip and Heatmap or on the right of Heatmap is ultra clean! Let's put gene names on ax_a_left with rdend hidden/placed cleanly, OR place gene labels cleanly on the right of ax_a_mat / left of ax_a_rbar?
# Actually in gs_a: [2,0]=rdend, [2,1]=Source Strip, [2,2]=Heatmap, [2,3]=Row Bar (with gene names on its left? No, space!).
# Even cleaner: put gene names directly inside ax_a_rbar's y-ticks on the RIGHT side (tick_right)!! Then rdend connects straight into Source Strip -> Heatmap -> Row Bar -> Italic Gene Name!! ZERO text colliding with dendrogram!
src_map = dict(zip(prot_df["gene"], prot_df["source"]))
for r_i, g in enumerate(gene_order):
    is_xray = str(src_map.get(g, "")).startswith("RCSB")
    col_s = NAVY if is_xray else "#8FA8BE"
    ax_a_left.add_patch(patches.Rectangle((0, r_i - 0.5), 1.0, 1.0, fc=col_s, ec="white", lw=1.0))
ax_a_left.set_xlim(0, 1); ax_a_left.set_ylim(len(gene_order) - 0.5, -0.5)
ax_a_left.axis("off")

# Luxurious Sequential Ivory -> Sky Teal -> Deep Sapphire (strong negative = rich dark navy-sapphire!)
vina_cmap = LinearSegmentedColormap.from_list(
    "sapphire_bind", ["#0B2545", "#1D5B9B", "#5BA4CF", "#C7DFEE", "#FAF8F5"]
)
im = ax_a_mat.imshow(mat_piv.values, cmap=vina_cmap, vmin=-10.2, vmax=-2.0, aspect="auto")

ax_a_mat.set_xticks(np.arange(-0.5, len(met_order), 1), minor=True)
ax_a_mat.set_yticks(np.arange(-0.5, len(gene_order), 1), minor=True)
ax_a_mat.grid(which="minor", color="white", linewidth=1.3)
ax_a_mat.tick_params(which="minor", bottom=False, left=False)

for r in range(len(gene_order)):
    for c in range(len(met_order)):
        val = mat_piv.values[r, c]
        ev  = str(ev_piv.values[r, c]).strip()
        has_ev = (ev != "")
        tcol = "white" if val <= -6.8 else INK
        txt = f"{val:.1f}*" if has_ev else f"{val:.1f}"
        ax_a_mat.text(c, r, txt, ha="center", va="center",
                      fontsize=7.0, fontweight="bold" if (has_ev or val <= -6.5) else "normal", color=tcol)
        if has_ev:
            ax_a_mat.add_patch(patches.Rectangle((c - 0.44, r - 0.44), 0.88, 0.88,
                                                 fill=False, ec="#F59E0B", lw=1.6, zorder=4))

ax_a_mat.set_xticks(range(len(met_order)))
ax_a_mat.set_xticklabels([met_disp[m] for m in met_order], rotation=32, ha="right", fontsize=7.6, color=INK)
ax_a_mat.set_yticks([])
for sp in ax_a_mat.spines.values(): sp.set_visible(False)

# Right marginal barplot: Strong pockets (<= -6.5 kcal/mol) + gene labels on right
strong_cnts = [(mat_piv.loc[g].values <= -6.5).sum() for g in gene_order]
y_pos = np.arange(len(gene_order))
ax_a_rbar.barh(y_pos, strong_cnts, height=0.62, color=TEAL, edgecolor="white", lw=0.6)
for yi, cnt in zip(y_pos, strong_cnts):
    if cnt > 0:
        ax_a_rbar.text(cnt + 0.18, yi, str(cnt), va="center", ha="left", fontsize=6.8, fontweight="bold", color=NAVY)
ax_a_rbar.set_ylim(len(gene_order) - 0.5, -0.5)
ax_a_rbar.set_xlim(0, 8.2)
ax_a_rbar.set_xticks([0, 4, 8])
ax_a_rbar.set_xlabel("# \u2264\u22126.5", fontsize=7.0, fontweight="bold", color=INK)
ax_a_rbar.yaxis.tick_right()
ax_a_rbar.set_yticks(y_pos)
ax_a_rbar.set_yticklabels(gene_order, fontsize=7.9, fontweight="bold", fontstyle="italic", color=INK)
ax_a_rbar.tick_params(axis="y", length=0, pad=3)
ax_a_rbar.spines["top"].set_visible(False)
ax_a_rbar.spines["left"].set_visible(False)
ax_a_rbar.spines["right"].set_visible(False)

# Horizontal colorbar cleanly below Panel A's x-labels + legend
cb_ax = ax_a_mat.inset_axes([0.0, -0.31, 0.39, 0.045])
cb6a = fig.colorbar(im, cax=cb_ax, orientation="horizontal")
cb6a.set_label("Vina Affinity (kcal/mol)", fontsize=7.1, fontweight="bold", color=INK, labelpad=2)
cb6a.set_ticks([-10, -8, -6.5, -4, -2])
cb6a.ax.tick_params(labelsize=6.5, pad=1)
cb6a.outline.set_linewidth(0.45)

leg_handles_a = [
    patches.Patch(fc=NAVY, label="RCSB X-ray"),
    patches.Patch(fc="#8FA8BE", label="AlphaFold"),
    patches.Patch(fc="none", ec="#F59E0B", lw=1.4, label="*gutMGene2"),
    patches.Patch(fc="#8C6D46", label="Bile Acid"),
    patches.Patch(fc="#1F8A70", label="Polyphenol"),
    patches.Patch(fc="#3C5488", label="Indole"),
    patches.Patch(fc="#2B7A9E", label="SCFA"),
    patches.Patch(fc="#7E8B9B", label="B3 / Choline"),
]
ax_a_mat.legend(handles=leg_handles_a, loc="upper left", bbox_to_anchor=(0.43, -0.22),
                ncol=3, frameon=False, fontsize=6.6, handlelength=1.0, columnspacing=0.8)

# ==============================================================================
# PANEL B: Binding Energy Spectrum by Scaffold (Raincloud + Exact Callouts)
# ==============================================================================
ax6b = fig.add_subplot(top_gs[1])
class_map_b = {
    "Lithocholic_acid": "Bile Acid",
    "Urolithin_A": "Polyphenol /\nPhenolic",
    "Caffeic_acid": "Polyphenol /\nPhenolic",
    "3-(4-Hydroxyphenyl)propionic_acid": "Polyphenol /\nPhenolic",
    "3-Indolepropionic_acid": "Tryptophan\nIndoles",
    "Indole-3-carboxylic_acid": "Tryptophan\nIndoles",
    "Indole": "Tryptophan\nIndoles",
    "Nicotinic_acid": "Pyridine /\nB3",
    "Butyrate": "SCFA",
    "Propionate": "SCFA",
    "Trimethylamine_N-oxide": "TMAO",
}
dock_df["Scaffold"] = dock_df["metabolite"].map(class_map_b)
scaf_ord  = ["Bile Acid", "Polyphenol /\nPhenolic", "Tryptophan\nIndoles", "Pyridine /\nB3", "SCFA", "TMAO"]
scaf_cols = ["#8C6D46", "#1F8A70", "#3C5488", "#7E8B9B", "#2B7A9E", "#C84B31"]

rng = np.random.default_rng(2026)
pt_coords = {}
for idx, sc in enumerate(scaf_ord):
    sub_df = dock_df[dock_df["Scaffold"] == sc].reset_index(drop=True)
    vals = sub_df["affinity_kcal_mol"].values
    col  = scaf_cols[idx]
    parts = ax6b.violinplot([vals], positions=[idx - 0.05], widths=0.55, showmeans=False, showmedians=False, showextrema=False)
    for pc in parts["bodies"]:
        m = np.mean(pc.get_paths()[0].vertices[:, 0])
        pc.get_paths()[0].vertices[:, 0] = np.clip(pc.get_paths()[0].vertices[:, 0], -np.inf, m)
        pc.set_facecolor(col); pc.set_alpha(0.30); pc.set_edgecolor(col); pc.set_linewidth(0.75)
    ax6b.boxplot([vals], positions=[idx + 0.06], widths=0.18, patch_artist=True, showfliers=False,
                 boxprops=dict(facecolor=col, alpha=0.85, edgecolor=INK, lw=0.7),
                 medianprops=dict(color="white", lw=1.3),
                 whiskerprops=dict(color=INK, lw=0.7), capprops=dict(color=INK, lw=0.7))
    jxs = idx + 0.22 + rng.uniform(-0.05, 0.05, len(vals))
    ax6b.scatter(jxs, vals, s=16, color=col, edgecolor="white", lw=0.45, alpha=0.90, zorder=4)
    for k_i, row_k in sub_df.iterrows():
        pt_coords[(row_k["gene"], row_k["metabolite"])] = (jxs[k_i], row_k["affinity_kcal_mol"])

ax6b.axhline(-6.5, ls="--", color=CORAL, lw=0.95, zorder=1)
ax6b.text(5.48, -6.38, "Strong pocket threshold (\u22126.5 kcal/mol)", ha="right", va="bottom",
          fontsize=7.2, fontweight="bold", color=CORAL)

callouts = [
    (("MUC1", "Lithocholic_acid"),        "MUC1\u2013LCA (\u221210.05)", (0.95, -10.15)),
    (("PEMT", "Urolithin_A"),             "PEMT\u2013UroA (\u22127.84)",  (2.22, -9.05)),
    (("IDO1", "Indole-3-carboxylic_acid"), "IDO1\u2013I3CA (\u22126.68)", (3.18, -7.75)),
]
for key, txt, (tx, ty) in callouts:
    px, py = pt_coords[key]
    ax6b.annotate(txt, xy=(px, py), xytext=(tx, ty), fontsize=6.9, fontweight="bold", color=INK,
                  arrowprops=dict(arrowstyle="-", color="#4A5568", lw=0.7),
                  bbox=dict(boxstyle="round,pad=0.18", fc="#F8FAFC", ec="#CBD5E0", lw=0.6))

ax6b.set_xticks(range(len(scaf_ord)))
ax6b.set_xticklabels(scaf_ord, fontsize=7.6, fontweight="bold", color=INK)
ax6b.set_ylabel("Vina Binding Energy (kcal/mol)", fontweight="bold", color=INK)
ax6b.set_ylim(-10.8, -1.5)
ax6b.set_xlim(-0.55, 5.65)
ax6b.spines["top"].set_visible(False); ax6b.spines["right"].set_visible(False)
add_panel_header(ax6b, "B", "Binding Energy by Metabolite Scaffold", y_axes=1.06, x_letter_pts=-18)

# ==============================================================================
# PANEL C: Borderless Supersampled 3D Pockets + Integrated Chemical Structures
# ==============================================================================
key_cards = [
    ("PEMT",  "Urolithin_A",              "PEMT \u00D7 Urolithin A",     "AlphaFold | \u22127.84 kcal/mol", "C1=CC2=C(C=C1O)C(=O)OC3=C2C=CC(=C3)O", "Urolithin A  (C13H8O4)"),
    ("IDO1",  "Indole-3-carboxylic_acid", "IDO1 \u00D7 Indole-3-CA",     "AlphaFold | \u22126.68 kcal/mol", "C1=CC=C2C(=C1)C(=CN2)C(=O)O",          "Indole-3-CA  (C9H7NO2)"),
    ("DAPK3", "Urolithin_A",              "DAPK3 \u00D7 Urolithin A",    "PDB 3BHY (1.24 \u00C5) | \u22127.65 kcal/mol", "C1=CC2=C(C=C1O)C(=O)OC3=C2C=CC(=C3)O", "Urolithin A  (C13H8O4)"),
    ("CASP1", "Caffeic_acid",             "CASP1 \u00D7 Caffeic acid",   "PDB 8WRA (1.45 \u00C5) | \u22126.59 kcal/mol", "C1=CC(=C(C=C1/C=C/C(=O)O)O)O",         "Caffeic acid  (C9H8O4)"),
]

for p_i, (gene, met, title_s, sub_s, smi, lig_lbl) in enumerate(key_cards):
    cell_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=mid_gs[p_i], height_ratios=[0.85, 0.31], hspace=0.10)
    ax_3d = fig.add_subplot(cell_gs[0])
    ax_2d = fig.add_subplot(cell_gs[1])

    img_3d = Image.open(os.path.join(POCKET_DIR, f"{gene}_{met}.png"))
    ax_3d.imshow(img_3d)
    ax_3d.axis("off")
    ax_3d.set_title(f"{title_s}\n{sub_s}", fontsize=8.2, fontweight="bold", color=INK, pad=4)

    lig_2d = render_rdkit_png(smi, width=380, height=120)
    ax_2d.imshow(lig_2d)
    ax_2d.axis("off")
    # Subtle top hairline and clean title ABOVE molecule so zero atom collision
    ax_2d.plot([0.06, 0.94], [1.14, 1.14], transform=ax_2d.transAxes, color="#E2E8F0", lw=0.8, clip_on=False)
    ax_2d.set_title(lig_lbl, fontsize=6.9, fontweight="bold", color=MUTED, pad=4)

    if p_i == 0:
        add_panel_header(ax_3d, "C", "3D Binding Pocket Architecture & Polar H-Bonds (\u22643.55 \u00C5)",
                         y_axes=1.15, x_letter_pts=-6)

# ==============================================================================
# PANEL D: Editorial Biomedical Illustration + Crisp Vector Data Overlay
# ==============================================================================
ax_d = fig.add_subplot(outer[2])
bg_full = Image.open(BG_6D_PATH).convert("RGB")
# Crop outer top/bottom whitespace tightly while keeping all illustration & generous annotation margins
# Content spans x:20..1350, y:141..641
bg_crop = bg_full.crop((0, 95, 1376, 685))
ax_d.imshow(bg_crop, extent=[0, 100, 3.0, 29.5], aspect="auto", zorder=0)
ax_d.set_xlim(0, 100); ax_d.set_ylim(0, 35.5)
ax_d.axis("off")
add_panel_header(ax_d, "D", "Integrated Gut\u2013Host\u2013TNBC Translational Mechanism", y_axes=1.01, x_letter_pts=0)

# Zone headers cleanly above the 3 compartments
ax_d.text(17.5, 32.8, "1. Gut Dysbiosis & Metabolic Rewiring", ha="center", va="center",
          fontsize=8.6, fontweight="bold", color=TEAL)
ax_d.text(49.5, 32.8, "2. Systemic Flux, Causal MR & Docking", ha="center", va="center",
          fontsize=8.6, fontweight="bold", color=NAVY)
ax_d.text(82.5, 32.8, "3. TNBC Microenvironment & ICI Response", ha="center", va="center",
          fontsize=8.6, fontweight="bold", color=CORAL)

# Zone 1 overlays (Gut Lumen & Mucosa, x=0..34)
ax_d.text(8.8, 29.8, "Commensal SCFA Producers \u2193", ha="center", va="center", fontsize=6.9, fontweight="bold", color=TEAL)
ax_d.text(8.8, 27.6, "Lachnospira, Roseburia, Dorea\n(Pure TNBC MR OR = 0.75)",
          ha="center", va="center", fontsize=6.2, style="italic", color=INK)

ax_d.text(25.2, 27.9, "Pathobionts \u2191", ha="center", va="center", fontsize=6.9, fontweight="bold", color=CORAL)
ax_d.text(25.2, 25.8, "Eggerthella, Sellimonas\n(TNBC MR OR = 1.12)",
          ha="center", va="center", fontsize=6.2, style="italic", color=INK)

ax_d.text(17.2, 2.2, "8 HGMT Cohorts (N = 1,629 runs): Healthy \u2192 Benign \u2192 BC gradient\n"
                     "SCFAs & Urolithin A \u2193   |   WGS Steroid (+1.13) & Tryptophan (+0.26) \u2191",
          ha="center", va="center", fontsize=6.3, fontweight="bold", color=INK)

# Zone 2 overlays (Blood vessel x=35..54, 3 stacked receptors x=55..65)
ax_d.text(45.5, 22.8, "Portal / Systemic Metabolite Flux", ha="center", va="center",
          fontsize=7.0, fontweight="bold", color=CORAL)
ax_d.text(45.5, 20.9, "Urolithin A \u2022 Indole-3-CA \u2022 Caffeic acid \u2022 SCFAs", ha="center", va="center",
          fontsize=6.1, color=INK)

# Labels cleanly tucked beside the 3 receptors BEFORE the x=65.3 divider line
ax_d.text(59.6, 27.6, "PEMT (Protective)", ha="left", va="center", fontsize=6.4, fontweight="bold", color=EMER)
ax_d.text(59.6, 25.6, "TNBC FDR=0.048\nUroA \u22127.84 kcal/mol", ha="left", va="center", fontsize=5.7, color=INK)

ax_d.text(60.5, 18.6, "IDO1 (Immune)", ha="left", va="center", fontsize=6.4, fontweight="bold", color="#5B4B8A",
          bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.88))
ax_d.text(60.5, 17.0, "I3CA \u22126.68", ha="left", va="center", fontsize=5.7, color=INK,
          bbox=dict(boxstyle="round,pad=0.06", fc="white", ec="none", alpha=0.88))

ax_d.text(59.6, 10.6, "DAPK3 / AKT1", ha="left", va="center", fontsize=6.4, fontweight="bold", color=TEAL)
ax_d.text(59.6,  8.7, "UroA \u22127.65 kcal/mol\nAKT1 GTEx P=0.006", ha="left", va="center", fontsize=5.7, color=INK)

ax_d.text(46.5, 2.2, "Steiger Forward: 56/56 TRUE  |  Reverse MR Confounding: 0/56 (FDR < 0.05)",
          ha="center", va="center", fontsize=6.2, fontweight="bold", color=NAVY)

# Zone 3 Top: Low GMSig-4 (Immune-Hot) - placed cleanly in empty white space ABOVE cells!
ax_d.text(66.5, 30.3, "Low GMSig-4 | Immune-Hot Niche",
          ha="left", va="center", fontsize=7.0, fontweight="bold", color=EMER)
ax_d.text(98.5, 30.3, "Pembro pCR: 78.6% | AUC = 0.737", ha="right", va="center",
          fontsize=6.1, fontweight="bold", color=EMER,
          bbox=dict(boxstyle="round,pad=0.20", fc="#ECF7F4", ec=EMER, lw=0.8))
ax_d.text(72.2, 27.9, "DC / Myeloid (CXCL10 \u2191)", ha="center", va="center", fontsize=5.9, fontweight="bold", color=TEAL)
ax_d.text(86.2, 27.9, "CD8+ T (CXCR3 / IFNG)", ha="center", va="center", fontsize=5.9, fontweight="bold", color=EMER)

# Zone 3 Bottom: High GMSig-4 (Immune-Cold / Fibrotic Excluded)
ax_d.text(66.5, 15.8, "High GMSig-4 | Glycolytic & CAF-Excluded",
          ha="left", va="center", fontsize=7.0, fontweight="bold", color=CORAL)
ax_d.text(83.2, 4.1, "ENO1+ Tumor (Glycolysis +2.03, EMT +2.43)\niCAF CXCL12 & myCAF COL1A1 Barrier",
          ha="center", va="center", fontsize=5.8, fontweight="bold", color=CORAL)
ax_d.text(69.8, 4.1, "Excluded\nT Cells", ha="center", va="center", fontsize=5.8, fontweight="bold", color=MUTED)
ax_d.text(98.5, 4.1, "ICI Resistance\npCR: 53.3%", ha="right", va="center",
          fontsize=6.1, fontweight="bold", color=CORAL,
          bbox=dict(boxstyle="round,pad=0.20", fc="#FDEEEC", ec=CORAL, lw=0.8))

out_png = os.path.join(FIG_DIR, "Figure6_v2_pilot.png")
out_pdf = os.path.join(FIG_DIR, "Figure6_v2_pilot.pdf")
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print("SAVED flagship Figure6_v2_pilot.png and Figure6_v2_pilot.pdf successfully!")
