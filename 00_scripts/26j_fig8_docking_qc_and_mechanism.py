import os
import io
import pandas as pd
import numpy as np
import matplotlib
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

for fpath in [
    "/mnt/c/Windows/Fonts/arial.ttf",
    "/mnt/c/Windows/Fonts/arialbd.ttf",
    "/mnt/c/Windows/Fonts/ariali.ttf",
    "/mnt/c/Windows/Fonts/arialbi.ttf",
]:
    if os.path.exists(fpath):
        fm.fontManager.addfont(fpath)

plt.rcParams.update({
    "font.family": "Arial",
    "font.sans-serif": ["Arial", "DejaVu Sans"],
    "axes.edgecolor": "#2D3748",
    "axes.linewidth": 0.65,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.labelsize": 8.5,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
})

RES_DIR    = "/mnt/d/Bioinformatics/TNBC_Microbiome/02_results"
QC_DIR     = os.path.join(RES_DIR, "JTM_upgrade/docking")
FIG_DIR    = "/mnt/d/Bioinformatics/TNBC_Microbiome/03_figures"
POCKET_DIR = "/mnt/c/Users/Administrator/.gemini/antigravity/brain/cbddd897-d490-44ba-a15e-fff716fe00ba/scratch/pockets3d"
BG_6D_PATH = "/mnt/c/Users/Administrator/.gemini/antigravity/brain/cbddd897-d490-44ba-a15e-fff716fe00ba/panel6d_bg_1795175998396.jpg"
os.makedirs(FIG_DIR, exist_ok=True)

INK   = "#1A202C"
MUTED = "#4A5568"
TEAL  = "#2B7A9E"
EMER  = "#1F8A70"
NAVY  = "#1E3A5F"
CORAL = "#C84B31"
OCHRE = "#C07A28"

# Exclude mis-targeted Chain-A PDB entries (ENO1 5NI9=HLA-DRA, TCF7L2 1JDH=CTNNB1) per structural QC audit!
dock_all = pd.read_csv(os.path.join(RES_DIR, "Docking_Vina_Affinities_Matrix.csv"))
dock_df  = dock_all[~dock_all["gene"].isin(["ENO1", "TCF7L2"])].copy()
prot_df  = pd.read_csv(os.path.join(RES_DIR, "Docking_Protein_Manifest.csv"))
qc_df    = pd.read_csv(os.path.join(QC_DIR, "redocking_and_seed_qc.csv"))

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
                textcoords="offset points", fontsize=9.8, fontweight="bold", color=INK, va="bottom", ha="left")

fig = plt.figure(figsize=(15.6, 17.2), dpi=300, facecolor="white")
outer = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1.02, 0.92, 0.94], hspace=0.58)
top_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[1.26, 1.12], wspace=0.38)
mid_gs = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[1], wspace=0.14)

# ==============================================================================
# PANEL A: QC-Filtered 9x11 Clustered Vina Matrix (99 Valid Target-Ligand Pairs)
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

row_lnk = linkage(pdist(mat_raw.values, metric="euclidean"), method="ward")
col_lnk = linkage(pdist(mat_raw.values.T, metric="euclidean"), method="ward")
d_row = dendrogram(row_lnk, orientation="left", ax=ax_a_rdend, no_labels=True, color_threshold=0, above_threshold_color="#64748B")
d_col = dendrogram(col_lnk, orientation="top",  ax=ax_a_cdend, no_labels=True, color_threshold=0, above_threshold_color="#64748B")
for ax_d in (ax_a_rdend, ax_a_cdend):
    for col_line in ax_d.collections:
        col_line.set_linewidth(0.85)
    ax_d.axis("off")

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
class_colors = {"Bile Acid": "#8C6D46", "Polyphenol": "#1F8A70", "Indole": "#3C5488", "B3 / Choline": "#7E8B9B", "SCFA": "#2B7A9E"}
met_disp = {
    "Lithocholic_acid": "Lithocholic acid", "Urolithin_A": "Urolithin A",
    "3-Indolepropionic_acid": "IPA", "Caffeic_acid": "Caffeic acid",
    "Indole-3-carboxylic_acid": "Indole-3-CA", "3-(4-Hydroxyphenyl)propionic_acid": "4-HPPA",
    "Indole": "Indole", "Nicotinic_acid": "Nicotinic acid",
    "Butyrate": "Butyrate", "Propionate": "Propionate", "Trimethylamine_N-oxide": "TMAO",
}

for c_i, m in enumerate(met_order):
    ax_a_top.add_patch(patches.Rectangle((c_i - 0.5, 0), 1.0, 1.0, fc=class_colors[class_map[m]], ec="white", lw=1.0))
ax_a_top.set_xlim(-0.5, len(met_order) - 0.5); ax_a_top.set_ylim(0, 1); ax_a_top.axis("off")

add_panel_header(ax_a_cdend, "A", "QC-Filtered Vina Matrix (9 Audited Targets \u00D7 11 Metabolites; 5NI9/1JDH Excluded)",
                 y_axes=1.12, x_letter_pts=-68)

src_map = dict(zip(prot_df["gene"], prot_df["source"]))
for r_i, g in enumerate(gene_order):
    is_xray = str(src_map.get(g, "")).startswith("RCSB")
    ax_a_left.add_patch(patches.Rectangle((0, r_i - 0.5), 1.0, 1.0, fc=NAVY if is_xray else "#8FA8BE", ec="white", lw=1.0))
ax_a_left.set_xlim(0, 1); ax_a_left.set_ylim(len(gene_order) - 0.5, -0.5); ax_a_left.axis("off")

vina_cmap = LinearSegmentedColormap.from_list("sapphire_bind", ["#0B2545", "#1D5B9B", "#5BA4CF", "#C7DFEE", "#FAF8F5"])
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
        ax_a_mat.text(c, r, f"{val:.1f}*" if has_ev else f"{val:.1f}", ha="center", va="center",
                      fontsize=7.2, fontweight="bold" if (has_ev or val <= -6.5) else "normal", color=tcol)
        if has_ev:
            ax_a_mat.add_patch(patches.Rectangle((c - 0.44, r - 0.44), 0.88, 0.88, fill=False, ec="#F59E0B", lw=1.6, zorder=4))

ax_a_mat.set_xticks(range(len(met_order)))
ax_a_mat.set_xticklabels([met_disp[m] for m in met_order], rotation=32, ha="right", fontsize=7.6, color=INK)
ax_a_mat.set_yticks([])
for sp in ax_a_mat.spines.values(): sp.set_visible(False)

strong_cnts = [(mat_piv.loc[g].values <= -6.5).sum() for g in gene_order]
y_pos = np.arange(len(gene_order))
ax_a_rbar.barh(y_pos, strong_cnts, height=0.60, color=TEAL, edgecolor="white", lw=0.6)
for yi, cnt in zip(y_pos, strong_cnts):
    if cnt > 0:
        ax_a_rbar.text(cnt + 0.18, yi, str(cnt), va="center", ha="left", fontsize=6.9, fontweight="bold", color=NAVY)
ax_a_rbar.set_ylim(len(gene_order) - 0.5, -0.5); ax_a_rbar.set_xlim(0, 8.2); ax_a_rbar.set_xticks([0, 4, 8])
ax_a_rbar.set_xlabel("# \u2264\u22126.5", fontsize=7.0, fontweight="bold", color=INK)
ax_a_rbar.yaxis.tick_right()
ax_a_rbar.set_yticks(y_pos)
ax_a_rbar.set_yticklabels(gene_order, fontsize=8.0, fontweight="bold", fontstyle="italic", color=INK)
ax_a_rbar.tick_params(axis="y", length=0, pad=3)
for sp in ("top", "left", "right"): ax_a_rbar.spines[sp].set_visible(False)

cb_ax = ax_a_mat.inset_axes([0.0, -0.32, 0.39, 0.045])
cb6a = fig.colorbar(im, cax=cb_ax, orientation="horizontal")
cb6a.set_label("Vina Affinity (kcal/mol)", fontsize=7.1, fontweight="bold", color=INK, labelpad=2)
cb6a.set_ticks([-10, -8, -6.5, -4, -2]); cb6a.ax.tick_params(labelsize=6.5, pad=1); cb6a.outline.set_linewidth(0.45)

leg_handles_a = [
    patches.Patch(fc=NAVY, label="RCSB X-ray"), patches.Patch(fc="#8FA8BE", label="AlphaFold"),
    patches.Patch(fc="none", ec="#F59E0B", lw=1.4, label="*gutMGene2"),
    patches.Patch(fc="#8C6D46", label="Bile Acid"), patches.Patch(fc="#1F8A70", label="Polyphenol"),
    patches.Patch(fc="#3C5488", label="Indole"), patches.Patch(fc="#2B7A9E", label="SCFA"),
    patches.Patch(fc="#7E8B9B", label="B3 / Choline"),
]
ax_a_mat.legend(handles=leg_handles_a, loc="upper left", bbox_to_anchor=(0.43, -0.22),
                ncol=3, frameon=False, fontsize=6.6, handlelength=1.0, columnspacing=0.8)

# ==============================================================================
# PANEL B: Structural QC — Redocking RMSD (9FOH+A1IEE) & 5-Seed Repeatability
# ==============================================================================
gs_b = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=top_gs[1], height_ratios=[1.05, 1.0], hspace=0.45)
ax_b1 = fig.add_subplot(gs_b[0])
ax_b2 = fig.add_subplot(gs_b[1])

# Top subpanel B1: Positive-Control Redocking RMSD (apo-IDO1 9FOH Chain A + A1IEE)
redock_cats = [
    ("Focused 26\u00C5\n(ETKDGv3, 5/5)", [0.633, 0.640, 0.642, 0.645, 0.652], EMER),
    ("Focused 26\u00C5\n(Xtal 3D, 4/5)",  [0.565, 0.590, 0.610, 0.625, 8.590], TEAL),
    ("Focused 22\u00C5\n(Xtal 3D, 1/1)",  [0.603],                             NAVY),
    ("Blind Core Box\n(156k \u00C5\u00B3, FAIL)", [9.679],                     CORAL),
]
for xi, (lbl, vals, col) in enumerate(redock_cats):
    med_v = float(np.median(vals))
    ax_b1.bar(xi, med_v, width=0.52, color=col, alpha=0.28, edgecolor=col, lw=1.2)
    rng_j = np.linspace(-0.10, 0.10, len(vals)) if len(vals) > 1 else [0.0]
    ax_b1.scatter([xi + dx for dx in rng_j], vals, s=34, color=col, edgecolor="white", lw=0.6, zorder=5)
    txt_y = max(vals) + 0.45
    ax_b1.text(xi, min(txt_y, 10.6), f"{med_v:.2f} \u00C5" if len(vals) == 1 else f"med {med_v:.2f} \u00C5",
               ha="center", va="bottom", fontsize=7.0, fontweight="bold", color=col)

ax_b1.axhline(2.0, ls="--", color=CORAL, lw=1.0)
ax_b1.text(3.42, 2.25, "2.0 \u00C5 redocking reference", ha="right", va="bottom", fontsize=6.9, fontweight="bold", color=CORAL)
ax_b1.set_xticks(range(len(redock_cats)))
ax_b1.set_xticklabels([c[0] for c in redock_cats], fontsize=7.2, fontweight="bold", color=INK)
ax_b1.set_ylabel("Sym-Corrected\nHeavy-Atom RMSD (\u00C5)", fontsize=7.6, fontweight="bold", color=INK)
ax_b1.set_ylim(0, 11.4); ax_b1.set_xlim(-0.55, 3.55)
ax_b1.spines["top"].set_visible(False); ax_b1.spines["right"].set_visible(False)
add_panel_header(ax_b1, "B", "Structural QC: Apo-IDO1 Redocking (9FOH) & Multi-Seed Stability", y_axes=1.10, x_letter_pts=-18)

# Bottom subpanel B2: 5-Seed Repeatability & AlphaFold Pocket pLDDT
qc_bars = [
    ("PEMT\u00D7UroA\n(Blind)",     -7.849, 0.009, "RMSD 0.08\u00C5 | pLDDT 93.0", EMER),
    ("PEMT\u00D7UroA\n(26\u00C5)",  -7.872, 0.022, "RMSD 0.05\u00C5 | 5/5 Mode 1", TEAL),
    ("IDO1\u00D7I3CA\n(Blind)",     -6.682, 0.025, "2 modes (4.3\u00C5) | pLDDT 97.2", OCHRE),
    ("IDO1\u00D7I3CA\n(26\u00C5)",  -6.701, 0.008, "Mode B (0.03\u00C5) | HEM clash", CORAL),
]
yb = np.arange(len(qc_bars))
for yi, (lbl, mu, sd, note, col) in enumerate(qc_bars):
    ax_b2.barh(yi, abs(mu), xerr=sd, height=0.58, color=col, alpha=0.82, edgecolor="white", lw=0.7, capsize=2.5)
    ax_b2.text(abs(mu) + 0.16, yi, f"{mu:.2f}\u00B1{sd:.3f}  ({note})", va="center", ha="left", fontsize=6.7, fontweight="bold", color=INK)

ax_b2.set_yticks(yb)
ax_b2.set_yticklabels([x[0] for x in qc_bars], fontsize=7.2, fontweight="bold", color=INK)
ax_b2.invert_yaxis()
ax_b2.set_xlim(0, 12.6)
ax_b2.set_xlabel("|Vina Affinity| (kcal/mol; 5-Seed Mean \u00B1 SD)", fontsize=7.6, fontweight="bold", color=INK)
ax_b2.spines["top"].set_visible(False); ax_b2.spines["right"].set_visible(False)

# ==============================================================================
# PANEL C: 3D Pocket Architecture (with Explicit Apo / pLDDT QC Subtitles)
# ==============================================================================
key_cards = [
    ("PEMT",  "Urolithin_A",              "PEMT \u00D7 Urolithin A (Exploratory)",     "AF Q9UBM1 (in vacuo) | pLDDT 93.0 | \u22127.84", "C1=CC2=C(C=C1O)C(=O)OC3=C2C=CC(=C3)O", "Urolithin A  (C13H8O4)"),
    ("IDO1",  "Indole-3-carboxylic_acid", "apo-IDO1 \u00D7 Indole-3-CA (Heme-Free)",   "AF P14902 (vacated HEM) | pLDDT 97.2 | \u22126.68", "C1=CC=C2C(=C1)C(=CN2)C(=O)O",          "Indole-3-CA  (C9H7NO2)"),
    ("DAPK3", "Urolithin_A",              "DAPK3 \u00D7 Urolithin A",                  "PDB 3BHY (1.24 \u00C5, 26\u00C5 box) | \u22127.65", "C1=CC2=C(C=C1O)C(=O)OC3=C2C=CC(=C3)O", "Urolithin A  (C13H8O4)"),
    ("CASP1", "Caffeic_acid",             "CASP1 \u00D7 Caffeic acid",                 "PDB 8WRA (1.45 \u00C5) | \u22126.59 kcal/mol",      "C1=CC(=C(C=C1/C=C/C(=O)O)O)O",         "Caffeic acid  (C9H8O4)"),
]

for p_i, (gene, met, title_s, sub_s, smi, lig_lbl) in enumerate(key_cards):
    cell_gs = gridspec.GridSpecFromSubplotSpec(2, 1, subplot_spec=mid_gs[p_i], height_ratios=[0.85, 0.31], hspace=0.10)
    ax_3d = fig.add_subplot(cell_gs[0])
    ax_2d = fig.add_subplot(cell_gs[1])

    img_3d = Image.open(os.path.join(POCKET_DIR, f"{gene}_{met}.png"))
    ax_3d.imshow(img_3d); ax_3d.axis("off")
    ax_3d.set_title(f"{title_s}\n{sub_s}", fontsize=7.8, fontweight="bold", color=INK, pad=4)

    lig_2d = render_rdkit_png(smi, width=380, height=120)
    ax_2d.imshow(lig_2d); ax_2d.axis("off")
    ax_2d.plot([0.06, 0.94], [1.14, 1.14], transform=ax_2d.transAxes, color="#E2E8F0", lw=0.8, clip_on=False)
    ax_2d.set_title(lig_lbl, fontsize=6.9, fontweight="bold", color=MUTED, pad=4)

    if p_i == 0:
        add_panel_header(ax_3d, "C", "3D Pocket Hypotheses & Polar Contacts (\u22643.55 \u00C5; Rigid-Receptor In Vacuo)",
                         y_axes=1.15, x_letter_pts=-6)

# ==============================================================================
# PANEL D: Evidence-Graded Hypothesis-Generating Translational Framework
# ==============================================================================
ax_d = fig.add_subplot(outer[2])
bg_full = Image.open(BG_6D_PATH).convert("RGB")
bg_crop = bg_full.crop((0, 95, 1376, 685))
ax_d.imshow(bg_crop, extent=[0, 100, 3.0, 29.5], aspect="auto", zorder=0)
ax_d.set_xlim(0, 100); ax_d.set_ylim(0, 36.2); ax_d.axis("off")
add_panel_header(ax_d, "D", "Evidence-Graded Translational Framework (Solid: Replicated Cohorts | Dashed/Note: Exploratory Hypotheses)", y_axes=1.01, x_letter_pts=0)

ax_d.text(17.5, 33.4, "1. Upstream Microbiome Prior (Unpaired Cohorts)", ha="center", va="center", fontsize=8.3, fontweight="bold", color=TEAL)
ax_d.text(49.5, 33.4, "2. Genetic Prioritization & Structural Hypotheses", ha="center", va="center", fontsize=8.3, fontweight="bold", color=NAVY)
ax_d.text(82.5, 33.4, "3. Independently Validated TNBC Prognosis (GMSig-4)", ha="center", va="center", fontsize=8.3, fontweight="bold", color=CORAL)

ax_d.text(8.8, 29.8, "SCFA Producers \u2193 (Observed)", ha="center", va="center", fontsize=6.7, fontweight="bold", color=TEAL)
ax_d.text(8.8, 27.6, "Lachnospira, Roseburia, Dorea\n(Nominal TNBC MR OR = 0.75)", ha="center", va="center", fontsize=6.1, style="italic", color=INK)
ax_d.text(25.2, 27.9, "Pathobionts \u2191 (Observed)", ha="center", va="center", fontsize=6.7, fontweight="bold", color=CORAL)
ax_d.text(25.2, 25.8, "Eggerthella, Sellimonas\n(Nominal TNBC MR OR = 1.12)", ha="center", va="center", fontsize=6.1, style="italic", color=INK)
ax_d.text(17.2, 1.9, "8 HGMT Cohorts (1,629 runs): Healthy \u2192 Benign \u2192 BC gradient\nWGS Steroid (+1.13) & Tryptophan (+0.26) \u2191  |  Not patient-paired to tumor RNA", ha="center", va="center", fontsize=6.1, fontweight="bold", color=INK)

ax_d.text(45.5, 22.8, "Candidate Metabolite Hypotheses", ha="center", va="center", fontsize=6.8, fontweight="bold", color=CORAL)
ax_d.text(45.5, 20.9, "Urolithin A \u2022 Indole-3-CA \u2022 Caffeic acid \u2022 SCFAs", ha="center", va="center", fontsize=6.0, color=INK)

ax_d.text(59.6, 27.6, "PEMT (Blood MR)", ha="left", va="center", fontsize=6.3, fontweight="bold", color=EMER)
ax_d.text(59.6, 25.5, "TNBC FDR=0.048\nUroA \u22127.84 (in vacuo)", ha="left", va="center", fontsize=5.6, color=INK)
ax_d.text(60.5, 18.6, "apo-IDO1 (Hypoth.)", ha="left", va="center", fontsize=6.2, fontweight="bold", color="#5B4B8A",
          bbox=dict(boxstyle="round,pad=0.08", fc="white", ec="none", alpha=0.88))
ax_d.text(60.5, 16.9, "I3CA \u22126.68 (apo)", ha="left", va="center", fontsize=5.6, color=INK,
          bbox=dict(boxstyle="round,pad=0.06", fc="white", ec="none", alpha=0.88))
ax_d.text(59.6, 10.6, "DAPK3 / AKT1", ha="left", va="center", fontsize=6.3, fontweight="bold", color=TEAL)
ax_d.text(59.6,  8.6, "UroA \u22127.65 (3BHY)\nAKT1 GTEx P=0.006", ha="left", va="center", fontsize=5.6, color=INK)
ax_d.text(46.5, 1.9, "Steiger Fwd: 56/56 TRUE | Overlapping BCAC GWAS | Mammary coloc PP.H4 \u2264 0.11", ha="center", va="center", fontsize=6.0, fontweight="bold", color=NAVY)

ax_d.text(66.5, 30.3, "Low GMSig-4 | Immune-Hot", ha="left", va="center", fontsize=6.6, fontweight="bold", color=EMER)
ax_d.text(98.5, 30.3, "METABRIC OS P=0.0019 | GSE21653 P=6.6e-5", ha="right", va="center", fontsize=5.7, fontweight="bold", color=EMER,
          bbox=dict(boxstyle="round,pad=0.18", fc="#ECF7F4", ec=EMER, lw=0.8))
ax_d.text(72.2, 27.9, "CPTAC Immune Proteins \u2191", ha="center", va="center", fontsize=5.8, fontweight="bold", color=TEAL)
ax_d.text(86.2, 27.9, "STAT1/CD8A/GZMA (r=\u22120.86)", ha="center", va="center", fontsize=5.8, fontweight="bold", color=EMER)

ax_d.text(66.5, 15.8, "High GMSig-4 | Immune-Cold & Glycolytic (Poor OS)", ha="left", va="center", fontsize=6.8, fontweight="bold", color=CORAL)
gsea_hm = pd.read_csv(os.path.join(RES_DIR, "GMSig4_TCGA_TNBC_GSEA_Hallmarks.csv")).set_index("ID")["NES"]
gly_nes = float(gsea_hm.get("HALLMARK_GLYCOLYSIS", 2.00))
emt_nes = float(gsea_hm.get("HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION", 2.37))
ax_d.text(83.2, 4.1, f"ENO1+ Glycolysis (+{gly_nes:.2f}), EMT (+{emt_nes:.2f})\nExploratory ICI: Pembro 78.6% vs 53.3% (int P=0.14)", ha="center", va="center", fontsize=5.7, fontweight="bold", color=CORAL)
ax_d.text(69.8, 4.1, "Immune-Cold\nTME", ha="center", va="center", fontsize=5.8, fontweight="bold", color=MUTED)
ax_d.text(98.5, 4.1, "METABRIC +Clinical\nLRT P = 0.0089", ha="right", va="center", fontsize=5.9, fontweight="bold", color=CORAL,
          bbox=dict(boxstyle="round,pad=0.20", fc="#FDEEEC", ec=CORAL, lw=0.8))

out_png = os.path.join(FIG_DIR, "Figure8_Docking_QC_Mechanism.png")
out_pdf = os.path.join(FIG_DIR, "Figure8_Docking_QC_Mechanism.pdf")
fig.savefig(out_png, dpi=300, bbox_inches="tight", facecolor="white")
fig.savefig(out_pdf, bbox_inches="tight", facecolor="white")
print("SAVED upgraded Figure8_Docking_QC_Mechanism.png and .pdf successfully!")
