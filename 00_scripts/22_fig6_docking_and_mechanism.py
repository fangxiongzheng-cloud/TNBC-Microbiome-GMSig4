import os
import re
import numpy as np
import pandas as pd
import matplotlib

project_root = os.path.abspath(os.environ.get("PROJECT_ROOT") or ".")
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from matplotlib.patches import FancyBboxPatch

DATA_DIR = os.path.join(project_root, "01_data")
RES_DIR  = os.path.join(project_root, "02_results")
FIG_DIR  = os.path.join(project_root, "03_figures")
REC_DIR  = os.path.join(DATA_DIR, "docking", "prepared", "receptors")
POSE_DIR = os.path.join(DATA_DIR, "docking", "poses")
os.makedirs(FIG_DIR, exist_ok=True)

plt.rcParams.update({
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial"],
    "axes.edgecolor": "#444444",
    "axes.linewidth": 0.7,
    "xtick.labelsize": 8,
    "ytick.labelsize": 8,
    "axes.titlesize": 9.8,
    "axes.labelsize": 8.5,
})

# ------------------------------------------------------------
# Parsers for receptor PDB & AutoDock Vina MODEL 1 pose PDBQT
# ------------------------------------------------------------
def parse_receptor_pdb(path):
    atoms = []
    with open(path, "r") as f:
        for line in f:
            if line.startswith(("ATOM  ", "HETATM")):
                aname = line[12:16].strip()
                resn  = line[17:20].strip()
                chain = line[21:22].strip() or "A"
                resi  = line[22:26].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                elem = (line[76:78].strip() if len(line) >= 78 else "") or aname[0]
                if elem.upper().startswith("H"):
                    continue
                atoms.append((aname, resn, chain, resi, x, y, z, elem.upper()))
    return pd.DataFrame(atoms, columns=["aname", "resn", "chain", "resi", "x", "y", "z", "elem"])

def parse_ligand_model1(path):
    atoms = []
    in_m1 = False
    with open(path, "r") as f:
        for line in f:
            if line.startswith("MODEL 1"):
                in_m1 = True
                continue
            if in_m1 and line.startswith("ENDMDL"):
                break
            if in_m1 and line.startswith(("ATOM  ", "HETATM")):
                aname = line[12:16].strip()
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                atype = line[77:79].strip() if len(line) >= 79 else aname[0]
                elem = re.sub(r"[^A-Za-z]", "", atype)[0].upper() if atype else aname[0].upper()
                if elem == "H":
                    continue
                atoms.append((aname, x, y, z, elem, atype))
    return pd.DataFrame(atoms, columns=["aname", "x", "y", "z", "elem", "atype"])

def analyze_complex(gene, metabolite, cutoff=4.0):
    rec_df = parse_receptor_pdb(os.path.join(REC_DIR, f"{gene}_clean.pdb"))
    lig_df = parse_ligand_model1(os.path.join(POSE_DIR, f"{gene}__{metabolite}_out.pdbqt"))
    R = rec_df[["x", "y", "z"]].values
    L = lig_df[["x", "y", "z"]].values
    dists = np.linalg.norm(R[:, None, :] - L[None, :, :], axis=2)  # n_rec x n_lig
    min_r = dists.min(axis=1)
    mask = min_r <= cutoff
    near = rec_df.loc[mask].copy()
    near["min_dist"] = min_r[mask]
    # Flag potential polar contacts (N/O/S receptor to N/O/OA/NA/SA ligand <= 3.5 A)
    lig_polar = lig_df["elem"].isin(["N", "O", "S"]).values
    if lig_polar.any():
        polar_d = dists[mask][:, lig_polar].min(axis=1)
        near["polar_dist"] = np.where(near["elem"].isin(["N", "O", "S"]), polar_d, np.nan)
    else:
        near["polar_dist"] = np.nan

    res_sum = (near.sort_values("min_dist")
               .groupby(["resn", "chain", "resi"], as_index=False)
               .agg(min_dist=("min_dist", "min"),
                    closest_atom=("aname", "first"),
                    min_polar_dist=("polar_dist", "min"))
               .sort_values("min_dist"))
    res_sum["gene"] = gene
    res_sum["metabolite"] = metabolite
    res_sum["res_label"] = res_sum["resn"] + res_sum["resi"].astype(str)
    return rec_df, lig_df, res_sum

# Key causal / mechanistic complexes to analyze structurally
dock_df = pd.read_csv(os.path.join(RES_DIR, "Docking_Vina_Affinities_Matrix.csv"))
key_pairs = [
    ("PEMT",   "Urolithin_A",              "PEMT (AlphaFold) x Urolithin A\ngutMGene2 Activ. | TNBC MR FDR=0.048"),
    ("IDO1",   "Indole-3-carboxylic_acid", "IDO1 (AlphaFold) x Indole-3-CA\ngutMGene2 Inhib. (PMID:37770127)"),
    ("DAPK3",  "Urolithin_A",              "DAPK3 (3BHY 1.24A) x Urolithin A\ngutMGene2 Activ. | Co-crystal pocket"),
    ("CASP1",  "Caffeic_acid",             "CASP1 (8WRA 1.45A) x Caffeic acid\ngutMGene2 Inhib. (PMID:33096252)"),
]

contact_tables = []
complex_geom = {}
for gene, met, lbl in key_pairs:
    rec_df, lig_df, res_sum = analyze_complex(gene, met, cutoff=4.0)
    aff = float(dock_df[(dock_df["gene"]==gene) & (dock_df["metabolite"]==met)]["affinity_kcal_mol"].iloc[0])
    res_sum["affinity_kcal_mol"] = aff
    contact_tables.append(res_sum)
    complex_geom[(gene, met)] = (rec_df, lig_df, res_sum, aff, lbl)

contacts_all = pd.concat(contact_tables, ignore_index=True)
contacts_all.to_csv(os.path.join(RES_DIR, "Docking_Key_Complexes_Contact_Residues.csv"), index=False)
print("Saved Docking_Key_Complexes_Contact_Residues.csv; counts per pair:")
print(contacts_all.groupby(["gene", "metabolite"]).agg(n_res=("res_label", "count"),
      top_res=("res_label", lambda s: ",".join(list(s)[:5])),
      min_d=("min_dist", "min")))

# ------------------------------------------------------------
# Build Figure 6 (3 rows with generous gap under rotated x-labels)
# ------------------------------------------------------------
fig = plt.figure(figsize=(17.6, 15.0), dpi=300)
outer = gridspec.GridSpec(3, 1, figure=fig, height_ratios=[1.06, 1.05, 0.72], hspace=0.52)

top_gs = gridspec.GridSpecFromSubplotSpec(1, 2, subplot_spec=outer[0], width_ratios=[1.35, 1.0], wspace=0.25)
mid_gs = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=outer[1], wspace=0.22)
ax6d   = fig.add_subplot(outer[2])

# --- Panel 6A: 11x11 AutoDock Vina Affinity Matrix ---
ax6a = fig.add_subplot(top_gs[0])
mat_piv = dock_df.pivot(index="gene", columns="metabolite", values="affinity_kcal_mol")
ev_piv  = dock_df.fillna({"gutMGene_evidence": ""}).pivot(index="gene", columns="metabolite", values="gutMGene_evidence")

gene_order = mat_piv.mean(axis=1).sort_values().index.tolist()
met_order  = mat_piv.mean(axis=0).sort_values().index.tolist()
mat_piv = mat_piv.loc[gene_order, met_order]
ev_piv  = ev_piv.loc[gene_order, met_order]

met_short = {
    "Lithocholic_acid": "Lithocholic\nacid (BA)",
    "Urolithin_A": "Urolithin A\n(Polyphenol)",
    "3-Indolepropionic_acid": "IPA\n(Indole)",
    "Caffeic_acid": "Caffeic\nacid",
    "Indole-3-carboxylic_acid": "Indole-3-CA\n(Indole)",
    "3-(4-Hydroxyphenyl)propionic_acid": "4-HPPA\n(Phenolic)",
    "Indole": "Indole\n(Tryptophan)",
    "Nicotinic_acid": "Nicotinic\nacid (B3)",
    "Butyrate": "Butyrate\n(SCFA C4)",
    "Propionate": "Propionate\n(SCFA C3)",
    "Trimethylamine_N-oxide": "TMAO\n(Choline)",
}

im = ax6a.imshow(mat_piv.values, cmap="RdYlBu_r", vmin=-10.2, vmax=-2.2, aspect="auto")
ax6a.set_xticks(range(len(met_order)))
ax6a.set_xticklabels([met_short.get(m, m) for m in met_order], rotation=30, ha="right", fontsize=7.5)
ax6a.set_yticks(range(len(gene_order)))
ax6a.set_yticklabels(gene_order, fontsize=8.2, fontweight="bold")

for r in range(len(gene_order)):
    for c in range(len(met_order)):
        val = mat_piv.values[r, c]
        ev  = str(ev_piv.values[r, c]).strip()
        star = "*" if ev != "" else ""
        tcol = "white" if (val <= -7.6 or val >= -3.1) else "#111111"
        ax6a.text(c, r, f"{val:.2f}{star}", ha="center", va="center",
                  fontsize=7.0, fontweight="bold" if star or val <= -7.5 else "normal", color=tcol)
        if star:
            rect = plt.Rectangle((c - 0.48, r - 0.48), 0.96, 0.96, fill=False, ec="black", lw=1.2)
            ax6a.add_patch(rect)

cb6a = plt.colorbar(im, ax=ax6a, fraction=0.032, pad=0.02)
cb6a.set_label("Vina Affinity (kcal/mol)", fontsize=8.0)
ax6a.set_title("A. AutoDock Vina 1.2.7 Structural Docking Matrix (11 Causal Proteins x 11 Metabolites)\n(*black boxes denote curated gutMGene2 metabolite-target experimental links)",
               loc="left", fontweight="bold", fontsize=9.2)

# --- Panel 6B: Scaffold Affinity Distribution across 11 Targets ---
ax6b = fig.add_subplot(top_gs[1])
class_map = {
    "Lithocholic_acid": "Secondary Bile Acid",
    "Urolithin_A": "Ellagitannin / Phenolic",
    "Caffeic_acid": "Ellagitannin / Phenolic",
    "3-(4-Hydroxyphenyl)propionic_acid": "Ellagitannin / Phenolic",
    "3-Indolepropionic_acid": "Tryptophan Indoles",
    "Indole-3-carboxylic_acid": "Tryptophan Indoles",
    "Indole": "Tryptophan Indoles",
    "Nicotinic_acid": "Pyridine / Vitamin",
    "Butyrate": "Short-Chain Fatty Acids",
    "Propionate": "Short-Chain Fatty Acids",
    "Trimethylamine_N-oxide": "TMAO",
}
dock_df["Ligand_Class"] = dock_df["metabolite"].map(class_map)
class_order = ["Secondary Bile Acid", "Ellagitannin / Phenolic", "Tryptophan Indoles",
               "Pyridine / Vitamin", "Short-Chain Fatty Acids", "TMAO"]
class_labels_disp = ["Secondary\nBile Acid", "Ellagitannin /\nPhenolic", "Tryptophan\nIndoles",
                     "Pyridine /\nVitamin B3", "Short-Chain\nFatty Acids", "TMAO\n(Choline)"]
class_cols  = ["#7E6148", "#00A087", "#3C5488", "#8491B4", "#4DBBD5", "#E64B35"]

for idx, cl in enumerate(class_order):
    vals = dock_df.loc[dock_df["Ligand_Class"] == cl, "affinity_kcal_mol"].values
    bp = ax6b.boxplot([vals], positions=[idx], widths=0.52, patch_artist=True, showfliers=False,
                      boxprops=dict(facecolor=class_cols[idx], alpha=0.72, lw=0.8),
                      medianprops=dict(color="black", lw=1.4),
                      whiskerprops=dict(color="#444"), capprops=dict(color="#444"))
    jitter = np.linspace(-0.12, 0.12, len(vals))
    ax6b.scatter(idx + jitter, np.sort(vals), s=15, color=class_cols[idx], edgecolor="#222", lw=0.4, zorder=3)

ax6b.axhline(-6.5, ls="--", color="#B40426", lw=0.9, label="Strong direct pocket threshold (-6.5 kcal/mol)")
ax6b.set_xticks(range(len(class_order)))
ax6b.set_xticklabels(class_labels_disp, rotation=40, ha="right", fontsize=7.6)
ax6b.set_ylabel("Vina Binding Energy (kcal/mol)")
ax6b.set_ylim(-10.8, -1.5)
ax6b.grid(axis="y", color="#eeeeee", lw=0.5)
ax6b.legend(loc="lower right", fontsize=7.6, frameon=False)
ax6b.set_title("B. Binding Energy Spectrum by Microbial Metabolite Scaffold", loc="left", fontweight="bold", fontsize=9.2)

# --- Panel 6C (4 subpanels across Row 2): 3D Orthogonal Pocket Contact Diagrams ---
elem_colors = {"C": "#333333", "N": "#2B5B9E", "O": "#D62728", "S": "#DAA520"}
for p_idx, (gene, met, lbl) in enumerate(key_pairs):
    ax_p = fig.add_subplot(mid_gs[p_idx])
    rec_df, lig_df, res_sum, aff, _ = complex_geom[(gene, met)]
    L = lig_df[["x", "y", "z"]].values
    l_ctr = L.mean(axis=0)
    Lc = L - l_ctr
    # Orthogonal PCA projection along ligand principal axes for clearest flat pocket view
    _, _, vt = np.linalg.svd(Lc, full_matrices=False)
    proj = vt[:2].T  # 3x2
    L2 = Lc @ proj

    # Draw ligand covalent bonds (inter-atomic distance <= 1.85 A)
    ldist = np.linalg.norm(L[:, None, :] - L[None, :, :], axis=2)
    for i in range(len(L)):
        for j in range(i + 1, len(L)):
            if ldist[i, j] <= 1.85:
                ax_p.plot([L2[i, 0], L2[j, 0]], [L2[i, 1], L2[j, 1]], color="#111111", lw=2.2, zorder=4)

    for i, row in lig_df.reset_index(drop=True).iterrows():
        c = elem_colors.get(row["elem"], "#555555")
        ax_p.scatter(L2[i, 0], L2[i, 1], s=42 if row["elem"] in ["N", "O"] else 28,
                     c=c, edgecolors="white", linewidths=0.5, zorder=5)

    # Project top contacting residues (<= 4.0 A) using their closest contact atom or CA
    top_res = res_sum.head(7).copy()
    R_near_list = []
    for _, rr in top_res.iterrows():
        sub_atom = rec_df[(rec_df["resn"]==rr["resn"]) & (rec_df["resi"]==rr["resi"]) & (rec_df["chain"]==rr["chain"])]
        dmat = np.linalg.norm(sub_atom[["x", "y", "z"]].values[:, None, :] - L[None, :, :], axis=2)
        min_idx = np.unravel_index(np.argmin(dmat), dmat.shape)
        r_pt = (sub_atom[["x", "y", "z"]].values[min_idx[0]] - l_ctr) @ proj
        l_pt = L2[min_idx[1]]
        dist_val = float(dmat[min_idx])
        is_polar = pd.notna(rr["min_polar_dist"]) and rr["min_polar_dist"] <= 3.55
        R_near_list.append((rr["res_label"], r_pt[0], r_pt[1], l_pt[0], l_pt[1], dist_val, is_polar))

    for r_lbl, rx, ry, lx, ly, dval, is_pol in R_near_list:
        line_col = "#00A087" if is_pol else "#888888"
        ls_style = "-" if is_pol else "--"
        ax_p.plot([rx, lx], [ry, ly], color=line_col, ls=ls_style, lw=1.1 if is_pol else 0.75, zorder=2)
        ax_p.scatter(rx, ry, s=46, c="#FFF5E6" if is_pol else "#EFEFEF",
                     edgecolors=line_col, linewidths=1.1, zorder=3)

    # Sort residues by angle around center and repel angularly to guarantee ZERO label overlap
    angs = [np.arctan2(y, x) for _, x, y, *_ in R_near_list]
    order_a = np.argsort(angs)
    sorted_ang = np.array([angs[k] for k in order_a])
    n_res_p = len(order_a)
    if n_res_p > 1:
        min_sep = 2.0 * np.pi / max(n_res_p + 0.5, 6.0)
        for _ in range(40):
            for ii in range(n_res_p):
                nxt = (ii + 1) % n_res_p
                diff = sorted_ang[nxt] - sorted_ang[ii]
                if nxt == 0:
                    diff += 2.0 * np.pi
                if diff < min_sep:
                    push = 0.5 * (min_sep - diff)
                    sorted_ang[ii] -= push
                    sorted_ang[nxt] += push
    placed_coords = {}
    max_r = max([np.hypot(x, y) for _, x, y, *_ in R_near_list] + [float(np.hypot(L2[:,0], L2[:,1]).max())])
    tx_list, ty_list = [], []
    for slot, k in enumerate(order_a):
        r_lbl, rx, ry, lx, ly, dval, is_pol = R_near_list[k]
        line_col = "#00A087" if is_pol else "#888888"
        a_out = sorted_ang[slot]
        r_box = max(np.hypot(rx, ry) + 1.85, max_r * 0.92 + 1.75)
        tx, ty = r_box * np.cos(a_out), r_box * np.sin(a_out)
        tx_list.append(tx); ty_list.append(ty)
        ax_p.plot([rx, tx * 0.88], [ry, ty * 0.88], color="#AAAAAA", lw=0.6, zorder=1)
        tag = f"{r_lbl}\n({dval:.1f}A)"
        ax_p.text(tx, ty, tag, fontsize=6.4, fontweight="bold" if is_pol else "normal",
                  color="#005F50" if is_pol else "#333333", ha="center", va="center",
                  bbox=dict(boxstyle="round,pad=0.14", fc="white", ec=line_col, alpha=0.92, lw=0.65), zorder=6)

    pad = 1.9
    all_x = list(L2[:, 0]) + [x[1] for x in R_near_list] + tx_list
    all_y = list(L2[:, 1]) + [x[2] for x in R_near_list] + ty_list
    ax_p.set_xlim(min(all_x) - pad, max(all_x) + pad)
    ax_p.set_ylim(min(all_y) - pad, max(all_y) + pad)
    ax_p.set_xticks([]); ax_p.set_yticks([])
    n_tot = len(res_sum)
    ax_p.set_title(f"{lbl}\nVina dG = {aff:.2f} kcal/mol | {n_tot} residues <=4.0A",
                   fontsize=7.9, fontweight="bold", pad=5)
    if p_idx == 0:
        ax_p.text(0.0, 1.24, "C. Atomic Binding Pockets & Polar Contact Residues (Solid green = polar <=3.55A; Dashed grey = hydrophobic/vdW <=4.0A)",
                  transform=ax_p.transAxes, fontsize=9.4, fontweight="bold", ha="left")

# --- Panel 6D: Verified Multi-Omics Mechanistic Summary Axis ---
ax6d.set_xlim(0, 100); ax6d.set_ylim(0, 10)
ax6d.axis("off")
ax6d.set_title("D. Integrated Gut Microbiome -> Host Causal Target -> Tumor Immune Exclusion & Immunotherapy Response Axis",
               loc="left", fontweight="bold", fontsize=9.5, pad=6)

boxes = [
    (1.0,  1.2, 18.5, 7.6, "#EAF6F9", "#0086A8",
     "1. Microbiome Dysbiosis\n(8 HGMT Cohorts, N=1,606)\n"
     "- Benign->Malignant gradient\n  Lachnospira/Roseburia depleted\n  Eggerthella/Escherichia enriched\n"
     "- WGS: Butanoate/Tryptophan/Steroid\n  KEGG alterations across CN/FR/US"),
    (21.5, 1.2, 21.0, 7.6, "#EEF2FA", "#3C5488",
     "2. Subtype-Specific Causal MR\n(6 BC GWAS + eQTLGen + GTEx)\n"
     "- Dorea OR=0.75 (PureTNBC p=0.0028)\n- Sellimonas OR=1.12 (TNBC p=0.014)\n"
     "- Host Targets: PEMT (FDR=0.048),\n  CCR6 (p=0.028), MUC1, FFAR2,\n  TCF7L2, AKT1 (GTEx p=0.0064)"),
    (44.5, 1.2, 16.5, 7.6, "#FFF2EC", "#D95F02",
     "3. Metabolite Docking\n(Vina 1.2.7, 11x11 Matrix)\n"
     "- Urolithin A x PEMT (-7.84)\n- Indole-3-CA x IDO1 (-6.68)\n- Urolithin A x DAPK3 (-7.65)\n- Caffeic acid x CASP1 (-6.59)\n- LCA x MUC1/PEMT (-10.05/-9.62)"),
    (63.0, 1.2, 18.0, 7.6, "#FDEEEC", "#E64B35",
     "4. scRNA Niche & TME Switch\n(42,512 Cells + 3 Cohorts n=521)\n"
     "- GMSig-4: IFNG(T), IDO1/CXCL10\n  (Myeloid), ENO1(Tumor Epi 90.5%)\n"
     "- High-Risk: CD8 T & SCFA sensing\n  lost, 14/14 checkpoints repressed,\n  EMT +2.43 / Glycolysis +2.03"),
    (83.0, 1.2, 16.0, 7.6, "#ECF7F2", "#00A087",
     "5. Prognosis & ICI pCR\n(TCGA/GEO/I-SPY2/Durva)\n"
     "- 6/6 survival endpoints p<0.05\n  Nomogram C=0.866 / 0.785\n"
     "- Pembro TNBC pCR: 78.6% vs 53%\n  AUC=0.737 (UniBeta 0.795)\n  vs Paclitaxel Ctrl AUC=0.572"),
]

for bx, by, bw, bh, fc, ec, txt in boxes:
    patch = FancyBboxPatch((bx, by), bw, bh, boxstyle="round,pad=0.35,rounding_size=0.6",
                           facecolor=fc, edgecolor=ec, linewidth=1.6)
    ax6d.add_patch(patch)
    lines = txt.split("\n")
    ax6d.text(bx + bw/2, by + bh - 0.85, lines[0], ha="center", va="top",
              fontsize=8.2, fontweight="bold", color=ec)
    ax6d.text(bx + 0.55, by + bh - 2.2, "\n".join(lines[1:]), ha="left", va="top",
              fontsize=6.8, color="#222222", linespacing=1.3)

for ax_x in [19.85, 42.85, 61.35, 81.35]:
    ax6d.annotate("", xy=(ax_x + 1.35, 5.0), xytext=(ax_x - 0.05, 5.0),
                  arrowprops=dict(arrowstyle="-|>", color="#333333", lw=1.8, mutation_scale=13))

fig.suptitle("Figure 6. Structural Microbial Metabolite Docking and Integrated Causal-to-Clinical Translational Axis",
             fontsize=12.5, fontweight="bold", x=0.04, ha="left", y=0.988)

fig.savefig(os.path.join(FIG_DIR, "Figure6.pdf"), bbox_inches="tight")
fig.savefig(os.path.join(FIG_DIR, "Figure6.png"), dpi=300, bbox_inches="tight")
print("SAVED Figure6.pdf and Figure6.png successfully!")
