# Computational Structural Quality Audit: Protein–Ligand Docking (Figure 6)

**Scope:** Non-destructive computational structural audit of `00_scripts/21a_fetch_structures_and_ligands.py`, `00_scripts/21b_autodock_vina_pipeline.py`, `02_results/Docking_Protein_Manifest.csv`, `02_results/Docking_Vina_Affinities_Matrix.csv`, `02_results/Docking_Key_Complexes_Contact_Residues.csv`, and `04_manuscript/Manuscript_JTM.md`, supported by multi-seed repeatability, box-sensitivity, crystal redocking, AlphaFold pLDDT/UniProt topology mapping, and heme-retention feasibility checks using existing WSL bioinformatics tools.

**Generated QC Deliverables & Scripts:**
- Scripts: `00_scripts/26d_audit_receptors_plddt_topology.py`, `00_scripts/26d_vina_qc_seeds_box_redock.py` (pre-existing `00_scripts/26d_docking_fetch_qc.py` left untouched)
- Outputs:
  - `02_results/JTM_upgrade/docking/receptor_plddt_and_topology.csv` (602 residues: PEMT 1–199, IDO1 1–403)
  - `02_results/JTM_upgrade/docking/redocking_and_seed_qc.csv` (48 QC runs/tests)
  - `02_results/JTM_upgrade/docking/docking_audit_summary.md` (this report)
  - `02_results/JTM_upgrade/docking/figure6_caption_and_limitations.md` (drop-in manuscript text)
- Intermediate QC summary JSONs: `01_data/JTM_upgrade/docking/qc_work/part1_audit_metrics.json`, `01_data/JTM_upgrade/docking/qc_work/part2_qc_summary.json`

---

## 1. Installed Software Inventory

| Environment | Executable / Interpreter | Available Docking / Structural Packages | Absent CLI / AD4 Tools |
|---|---|---|---|
| **Windows host** | `Python 3.13.14` (`conda.bat` present) | None (`vina`, `meeko`, `rdkit`, `Bio` all absent) | `vina`, `obabel`, `autogrid4`, `autodock4`, `prepare_receptor4.py` |
| **WSL** (`bioinfo-singlecell-python`) | `/root/miniconda3/envs/bioinfo-singlecell-python/bin/python` (`Python 3.12.13`) | `vina 1.2.7` (Python bindings), `meeko 0.8.0` (`mk_prepare_receptor.py`, `mk_prepare_ligand.py`, `mk_export.py`), `rdkit 2026.09.2`, `biopython 1.88` | `vina` CLI, `obabel`, `autogrid4`, `autodock4`, `prepare_receptor4.py` |

**Manuscript vs. Code Discrepancy (Ligand/Receptor Prep & Protonation):**
`04_manuscript/Manuscript_JTM.md` (line 77) describes *"121 complexes prepared via OpenBabel/Meeko PDBQT conversion"*. In reality, `obabel` is neither installed nor called by `21a`/`21b`: ligands were generated from **neutral** PubChem `IsomericSMILES` using **RDKit ETKDGv3** (`randomSeed=2026`) + MMFF94/UFF minimization and **Meeko `MoleculePreparation`** (e.g., `Indole-3-carboxylic_acid.pdbqt` stores neutral carboxylic acid `REMARK SMILES O=C(O)c1c[nH]c2ccccc12` rather than physiological pH 7.4 carboxylate), and receptors were converted via **`mk_prepare_receptor.py`** (`--default_altloc A --allow_bad_res`, with `-b` raw-PDB fallback).

---

## 2. Receptor Origin, Chain Selection, HETATM Stripping, and Search-Box Audit

### 2.1 Receptor Origin and Chain-Selection Audit Across All 11 Targets

In `21a_fetch_structures_and_ligands.py`, RCSB queries selected the lowest-resolution X-ray entry (`<= 2.5 A`) and requested `https://files.rcsb.org/download/{best_pdb}.pdb`. When `.pdb` download failed (`if chosen == af_path or not os.path.exists(chosen):` line 60), the script silently fell back to AlphaFold without updating `best_pdb`/`resolution` in `Docking_Protein_Manifest.csv`. In `21b_autodock_vina_pipeline.py`, `clean_and_box` kept **only one chain** (`chosen_chain = "A"`, else first `ATOM` chain), dropped non-`(" ", "A")` altlocs, stripped all `HETATM` lines, and filtered AlphaFold models (`>= 80` `CA` atoms at `pLDDT >= 70.0`) to `pLDDT >= 60.0`.

| Gene | UniProt | `best_pdb` (Res) | Receptor Actually Docked (`source`) | Chains Raw -> Kept | Res Raw -> Clean -> PDBQT | Non-HOH HETATMs Stripped (`HOH`) | Search Box (`Mode`, Center, Size, Volume) | Critical Audit Finding |
|---|---|---|---|---|---|---|---|---|
| **IDO1** | P14902 | `9FOH` (1.60 A) | **`AlphaFold:P14902`** | `['A'] -> A` | `403 -> 379 -> 379` | None (AF) | `FoldedCore_Box` `(-5.70, 2.02, -1.81)`, `60.0x51.2x60.0`, **`184,273.2 A^3`** | **`9FOH` is mmCIF-only** (`pdb_format_compatible: "N"`), causing silent fallback to apo AlphaFold model; 24 residues (`1–9, 363–376, 380`) dropped (`pLDDT < 60`). |
| **PEMT** | Q9UBM1 | None | **`AlphaFold:Q9UBM1`** | `['A'] -> A` | `199 -> 192 -> 192` | None (AF) | `FoldedCore_Box` `(-0.98, -0.07, 1.74)`, `55.9x40.1x52.5`, **`117,462.1 A^3`** | No experimental structure; integral 4-helix ER/MAM membrane protein docked in vacuo without lipid bilayer or SAM cofactor (`98–100, 180–181`). |
| **ENO1** | P06733 | `5NI9` (1.33 A) | `RCSB_PDB:5NI9` | `['A','B','C'] -> A` | `386 -> 179 -> 179` | `MRD:8, URE:8, MPD:16` (`312 HOH`) | `CoCrystal_Pocket(MRD)` `(35.30, 9.21, -9.26)`, `26^3`, `17,576.0 A^3` | **Mis-targeted chain!** `5NI9` Chain A is **HLA-DRA** (`DR alpha chain`), Chain B is `DRB1-4`, and Chain C is only a 15-aa `alpha-enolase (326–340)` peptide. Rows labeled `ENO1` actually docked into **HLA-DRA** around buffer precipitant `MRD` (2-methyl-2,4-pentanediol). |
| **TCF7L2** | Q9NQB0 | `1JDH` (1.90 A) | `RCSB_PDB:1JDH` | `['A','B'] -> A` | `546 -> 508 -> 479` | None (`393 HOH`) | `FoldedCore_Box` `(-2.73, 8.94, 40.12)`, `60.0x48.4x60.0`, `174,128.4 A^3` | **Mis-targeted chain!** `1JDH` Chain A is **beta-catenin (`CTNNB1` 135–663)**, whereas `HTCF-4` (`TCF7L2` 12–49) is Chain B. Rows labeled `TCF7L2` actually docked into beta-catenin. |
| **MUC1** | P15941 | `6BSC` (1.30 A) | `RCSB_PDB:6BSC` | `['A','B'] -> A` | `105 -> 57 -> 57` | None (`92 HOH`) | `FoldedCore_Box` `(-14.90, 10.98, 3.05)`, `32.8x38.4x28.0`, `35,296.7 A^3` | Autoproteolytically cleaved SEA heterodimer split by keeping only N-terminal Chain A (`57/105` residues). |
| **FFAR2** | O15552 | None | `AlphaFold:O15552` | `['A'] -> A` | `330 -> 293 -> 293` | None (AF) | `FoldedCore_Box` `(0.43, 1.65, 3.32)`, `60.0x56.7x60.0`, `204,238.8 A^3` | GPCR AlphaFold model docked blindly without lipid bilayer (`37` terminal/loop residues dropped). |
| **TGFBR1** | P36897 | `9J9D` (1.34 A) | **`AlphaFold:P36897`** | `['A'] -> A` | `503 -> 433 -> 433` | None (AF) | `FoldedCore_Box` `(-5.12, -8.10, -8.27)`, `60.0x60.0x60.0`, `216,000.0 A^3` | Silent fallback from mmCIF-only `9J9D` (`1.34 A`) to full-length AlphaFold model; `70` residues dropped (`pLDDT < 60`). |
| **CASP1** | P29466 | `8WRA` (1.45 A) | `RCSB_PDB:8WRA` | `['A'] -> A` | `247 -> 247 -> 247` | `EDO:12` (`300 HOH`) | `FoldedCore_Box` `(15.32, -1.99, -4.46)`, `55.7x43.0x60.0`, `143,979.3 A^3` | Single-chain construct (`130..404`); cryoprotectant `EDO` excluded as solvent; blind box (`5.33x` `27,000 A^3`). |
| **NLRP3** | Q96P20 | `3QF2` (1.70 A) | `RCSB_PDB:3QF2` | `['A','B'] -> A` | `196 -> 90 -> 90` | None (`247 HOH`) | `FoldedCore_Box` `(-22.31, 29.53, 27.84)`, `36.8x31.3x32.4`, `37,305.0 A^3` | PYD/DAPIN domain only (`residues 3–112`), lacking NACHT/LRR domains. |
| **DAPK3** | O43293 | `3BHY` (1.24 A) | `RCSB_PDB:3BHY` | `['A'] -> A` | `268 -> 268 -> 248` | `CL:1, 7CP:25` (`303 HOH`) | **`CoCrystal_Pocket(7CP)`** `(14.40, 21.37, 29.35)`, `26^3`, `17,576.0 A^3` | Proper focused co-crystal ATP-cleft box (one of only two `26x26x26 A` runs); `20` residues dropped during `mk_prepare_receptor.py`. |
| **MAGI3** | Q5TCQ9 | `3SOE` (1.60 A) | `RCSB_PDB:3SOE` | `['A'] -> A` | `101 -> 101 -> 100` | `EDO:4` (`157 HOH`) | `FoldedCore_Box` `(4.53, 4.95, 24.56)`, `32.2x32.4x47.1`, `49,106.3 A^3` | Isolated 3rd PDZ domain (`residues 601–691`); `1` residue dropped in PDBQT. |

All 11 prepared PDBQT files had `pdbqt_hetatm_count = 0`. Nine of 11 targets exceeded AutoDock Vina's `27,000 A^3` recommended grid volume (`35,296.7` to `216,000.0 A^3`).

### 2.2 HETATM Stripping & Heme Cofactor Omission in IDO1

- **Why IDO1 lacked heme:** `AlphaFold:P14902` contains polypeptide atoms only (no `HEM` prosthetic group). Moreover, `21b` strips all `HETATM` lines and blacklists `"FE"` inside `SOLVENT_IONS`.
- **Crystal superposition proves direct heme-volume occupation:**
  - Aligning holo-IDO1 crystal structures **`2D0T`** (2.30 A, `373` CA aligned, CA RMSD `1.316 A`, bound `HEM` + distal `PIM`) and **`4PK5`** (2.79 A, `373` CA aligned, CA RMSD `1.359 A`, bound `HEM` + `PKJ`) onto `IDO1_P14902_AF` places `21b`'s top-ranked `Indole-3-carboxylic_acid` pose (`-6.677 kcal/mol`) directly inside the **vacated heme pocket**:
    - Minimum heavy-atom distance to aligned `2D0T HEM`: **`0.949 A`** (`4.329 A` to `HEM FE`).
    - Minimum heavy-atom distance to aligned `4PK5 HEM`: **`1.010 A`** (compared with `4.195 A` to `PKJ` and `4.979 A` to `2D0T PIM`).
    - Key pocket contacts in `Docking_Key_Complexes_Contact_Residues.csv` directly overlap heme-coordinating/heme-packing residues: **`HIS346`** (UniProt proximal axial `heme b` ligand, `1.8328 A` from `2D0T HEM`; docked `I3CA` `3.4334 A` to `CE1` / `3.5437 A` to `ND1`, `4.317 A` to `NE2`), **`PHE270`** (`0.2898 A` from `2D0T HEM`), **`ALA264`** (`1.7329 A`), and **`PHE214`** (`2.0632 A`).
  - Aligning apo-IDO1 crystal structure **`9FOH`** (1.60 A, `354` CA aligned, CA RMSD `1.368 A`, heme-free monomer complexed with `A1IEE` / `iDeg-2`, `inter_mol_metalic_bond_count: 0`) shows `I3CA` Pose 1 overlaps `A1IEE` at **`0.324 A`** minimum heavy-atom distance and shares 9 pocket residues (`TYR126, SER167, VAL170, GLU171, PHE214, ALA264, SER267, PHE270, HIS346`).
  - Consequently, the `21b` IDO1–I3CA score (`-6.68 kcal/mol`) models binding into a **heme-depleted apo cavity** and is **sterically forbidden in holo (heme-bound) IDO1** (`< 1.01 A` heavy-atom clash with `HEM`).

---

## 3. AlphaFold pLDDT Confidence & UniProt Topology Mapping (PEMT & IDO1)

Per-residue metrics are saved in `02_results/JTM_upgrade/docking/receptor_plddt_and_topology.csv` (602 rows).

### 3.1 PEMT (`Q9UBM1`, 199 residues)
- **Overall pLDDT:** Mean **`93.61`**, SD **`10.29`**, median **`96.56`** (IQR `94.34–98.00`), range **`37.25–98.69`**. Distribution: `177` very high (`>=90`), `12` confident (`70–90`), `6` low (`50–70`), `4` very low (`<50`).
- **Residues dropped by `21b` (`pLDDT < 60`):** `7` terminal residues removed (`1–2` [Lumenal, `45.44–54.69`]; `195–199` [Cytoplasmic, `37.25–58.72`]), leaving residues `3–194` (`192` kept).
- **UniProt topology & pLDDT spans:**
  - `1–12` Lumenal: mean `76.30` (`45.44–93.50`, `10/12` kept)
  - `13–33` Intramembrane Helical: mean `97.18` (`95.56–98.44`, `21/21` kept)
  - `34–45` Lumenal: mean `95.59` (`93.75–96.69`, `12/12` kept)
  - `46–66` Transmembrane Helical: mean `96.94` (`94.94–98.31`, `21/21` kept)
  - `67–93` Cytoplasmic: mean `96.46` (`91.94–98.56`, `27/27` kept)
  - `94–114` Transmembrane Helical: mean `97.52` (`94.56–98.62`, `21/21` kept)
  - `115–157` Lumenal: mean `95.58` (`87.94–98.69`, `43/43` kept)
  - `158–178` Transmembrane Helical: mean `97.52` (`95.44–98.56`, `21/21` kept)
  - `179–199` Cytoplasmic: mean `79.97` (`37.25–94.88`, `16/21` kept)
  - Combined TM/intramembrane helices (`84` residues): pLDDT mean **`97.29`**, min **`94.56`**.
- **Binding-pocket pLDDT & topology (`PEMT × Urolithin A`, `-7.836 kcal/mol`):**
  - All `10` Mode 1 contact residues (`<= 4.0 A`) are high confidence: pLDDT mean **`93.03`**, range **`87.94–96.25`** (`8` very high `>=90`; `2` confident `70–90`: `ARG129 = 87.94`, `LYS126 = 89.44`).
  - Topology breakdown of Mode 1 contacts: `LEU111` (`3.3086 A`, `95.56`, TM Helical `94–114`); `GLY115` (`3.5927 A`, `94.31`), `ASP120` (`3.2161 A`, `95.81`), `ILE124` (`3.9256 A`, `93.75`), `LYS126` (`3.7874 A`, `89.44`), `ARG129` (`2.9908 A`, `87.94`), `VAL130` (`3.7384 A`, `90.94`), `PHE135` (`3.4568 A`, `96.25`) (UniProt Lumenal `115–157`); and `THR184` (`3.3631 A`, `93.94`), `TYR188` (`3.2180 A`, `92.31`) (UniProt Cytoplasmic `179–199`).
  - **Topology & membrane/cofactor caveats:** In AlphaFold `Q9UBM1`, segment `136–157` folds across the Z-axis (`PHE135 Z=+14.46` to `HIS157 Z=-8.57`) and TM helix `158–178` returns (`ALA158 Z=-10.63` to `LEU178 Z=+0.70`), juxtaposing residues `111–135` and `180–188` (`GLU180` SAM-binding residue is **`4.1019 A`** from Mode 1 and **`4.1161 A`** from Mode 2; mutagenesis of `GLY98`, `GLY100`, `GLU180`, `GLU181` impairs or abolishes S-adenosyl-L-methionine binding). However, `21b` docks PEMT **in vacuo without the ER/MAM lipid bilayer or SAM cofactor**, and Mode 2 (`-7.764 kcal/mol`, `delta = 0.072 kcal/mol`) jumps **`12.308 A`** into the buried transmembrane core (`TRP29` [Intramembrane], `ASN64`, `VAL104`, `TYR173` [TM], `ARG67`, `THR116`, `PHE117`, `LEU118`, `TYR143`, `SER146`).

### 3.2 IDO1 (`P14902`, 403 residues)
- **Overall pLDDT:** Mean **`93.08`**, SD **`13.00`**, median **`97.81`** (IQR `95.34–98.50`), range **`33.38–98.94`**. Distribution: `350` very high (`>=90`), `22` confident (`70–90`), `18` low (`50–70`), `13` very low (`<50`).
- **Topology & residues dropped by `21b`:** Annotated `Cytoplasmic` (globular cytosolic enzyme; 0 TM helices). `21b` stripped `24` residues (`pLDDT < 60`): N-terminal tail `1–9` (`33.38–56.62`) and internal UniProt `Region:Disordered` (`360–381`, mean pLDDT `52.70`, `7/22` kept), dropping `363–376` and `380` (`35.66–56.53`). Crystal structures `2D0T`, `4PK5`, and `9FOH` likewise lack density across `361–379`, meaning `21b` introduced artificial peptide chain breaks (`PRO362 -> LYS377` and `THR379 -> GLY381`) leaving a disconnected 3-residue island (`377–379`, pLDDT `67.94–68.44`).
- **Binding-pocket pLDDT (`IDO1 × Indole-3-carboxylic acid`, `-6.677 kcal/mol`):** All `11` contact residues (`TYR126, SER167, VAL170, GLU171, PHE214, ALA264, SER267, VAL269, PHE270, LEU342, HIS346`) are folded-core cytoplasmic residues with pLDDT mean **`97.22`**, range **`88.88–98.69`** (`10 >= 90`; `ALA264 = 88.88`; `HIS346 = 96.94`). Thus, local backbone confidence in AlphaFold `P14902` is high; the critical defect is **omission of `heme b`** rather than low pocket pLDDT.

---

## 4. Multi-Seed Repeatability, Box Sensitivity, `9FOH` Redocking & Heme Feasibility (`redocking_and_seed_qc.csv`)

### 4.1 Multi-Seed Repeatability and Box-Size/Center Sensitivity

| System | Test Regime | Box Volume (`A^3`) | Seeds (`n`) | Pose 1 Affinity Mean +- SD (`Range`) kcal/mol | Symmetry-Corrected Pose 1 RMSD | Key Finding |
|---|---|---|---|---|---|---|
| **PEMT + UroA** | `21b_FoldedCore_Box` (blind, exh=16) | `117,462.1` | 5 | **`-7.8492 +- 0.0090`** (`-7.864 to -7.841`) | Pairwise mean **`0.081 A`** (max `0.133 A`); vs `21b`: `0.037–0.099 A` | **100% (`5/5`)** converged to Mode 1 (`<= 2.0 A`); Mode 2 persists within `0.074–0.096 kcal/mol`. |
| **PEMT + UroA** | `Focused_Pose1_26A` (exh=16) | `17,576.0` | 5 | **`-7.8724 +- 0.0223`** (`-7.892 to -7.836`) | Pairwise mean **`0.052 A`** (max `0.087 A`); vs `21b`: `0.052–0.128 A` | **100%** converged cleanly to Mode 1. |
| **PEMT + UroA** | Sensitivity (`exh=32`, `+10 A` box, `20 A` box, `Targeted_SAM_site_26A`) | `8,000–206,030.1` | 4 single runs | `-7.826` to `-7.887` | vs `21b`: **`0.056–0.240 A`** | Completely invariant across box volume (`8,000` to `206,030 A^3`) and SAM-targeted centering. |
| **IDO1 + I3CA** | `21b_FoldedCore_Box` (blind, exh=16) | `184,273.2` | 5 | **`-6.6824 +- 0.0247`** (`-6.705 to -6.645`) | Pairwise mean **`1.775 A`** (max **`4.363 A`**); **`60%` (`6/10`) pairs `<= 2.0 A`** | Seeds `2026/101/777/1337` stayed in Mode A (`0.041–0.095 A` vs `21b`, `0.904–0.967 A` to `2D0T HEM`); **`seed=42` flipped to Mode B (`4.363 A` RMSD, `2.022 A` centroid shift, `-6.645 kcal/mol`)**. |
| **IDO1 + I3CA** | `Focused_Pose1_26A` (exh=16) | `17,576.0` | 5 | **`-6.7014 +- 0.0079`** (`-6.707 to -6.689`) | Pairwise mean **`0.034 A`** (max `0.049 A`); vs `21b`: **`4.336–4.352 A`** | **5/5 (`100%`) flipped uniformly to Mode B** (`1.260–1.290 A` to `2D0T HEM`). |
| **IDO1 + I3CA** | `Expanded_plus10A` (`299,816.3 A^3`), `Focused_20A` (`8,000 A^3`), `Targeted_2D0T_PIM_distal_26A` (`17,576 A^3`) | `8,000–299,816.3` | 3 single runs (+`exh=32`) | `-6.593` to `-6.713` (`exh=32`: `-6.705`) | `exh=32`: `0.095 A` (Mode A); `+10 A`, `20 A`, `PIM_distal`: **`4.308–4.346 A` (Mode B)** | **Two isoenergetic apo orientations (`~0.01–0.04 kcal/mol` gap):** Mode A (`TYR126.OH 2.99 A / SER267.OG 2.82 A`) vs. Mode B (`ARG343.NE 3.07 A / HIS346.ND1 3.10 A`). Both clash with `2D0T HEM` (`0.904–1.671 A`). |

### 4.2 Crystal Ligand Redocking: Apo-IDO1 `9FOH` (Chain A) + `A1IEE` (iDeg-2)

Using heavy-atom symmetry-corrected RMSD (`rdkit.Chem.rdMolAlign.CalcRMS`) against the `9FOH.cif` experimental pose:

| Preparation / Box Regime | Volume (`A^3`) | Seeds (`n`) | Affinity Mean +- SD (`Range`) kcal/mol | Pose 1 Sym-RMSD vs Crystal | Best Top-5 RMSD Mean | Success Rate (`Pose 1 <= 2.0 A`) |
|---|---|---|---|---|---|---|
| `Pocket_26A` + `A1IEE_xtal_3D` | `17,576.0` | 5 | **`-12.454 +- 1.753`** (`-13.258 to -9.319`) | **`2.189 +- 3.578 A`** (`0.565–0.625 A` on 4 seeds; `8.590 A` on `seed=777`) | `1.498 A` (`0.565–0.625 A` vs `5.133 A`) | **`80%` (`4/5`)** |
| `Pocket_22A` + `A1IEE_xtal_3D` | `10,648.0` | 1 (`2026`) | `-13.286` | **`0.603 A`** | `0.603 A` | `100%` (`1/1`) |
| `Pocket_26A` + `A1IEE_SMILES_MMFF` (`21b` ETKDGv3 protocol) | `17,576.0` | 5 | **`-13.041 +- 0.018`** (`-13.055 to -13.012`) | **`0.642 +- 0.007 A`** (`0.633–0.652 A`) | `0.642 A` | **`100%` (`5/5`)** |
| **`9FOH_FoldedCore_Box` (blind `21b` protocol)** + `A1IEE_xtal_3D` | **`156,247.2`** | 1 (`2026`) | **`-9.398`** | **`9.679 A`** (best-of-10: `8.344 A`, centroid shift `4.420 A`) | `8.445 A` | **`0%` (FAILED)** |

**Interpretation:** Vina 1.2.7 + Meeko 0.8.0 (including `21b`'s exact `SMILES -> ETKDGv3/MMFF` ligand pipeline) accurately reproduces experimental IDO1 (`9FOH`) binding at **`0.63–0.65 A` RMSD** (`-13.04 kcal/mol`) when properly directed into the `26 A` pocket (though `1/5` `xtal_3D` seeds trapped in a `-9.32 kcal/mol` decoy), whereas `21b`'s automated `FoldedCore_Box` blind search (`156,247.2 A^3`) **completely failed** (`9.679 A` RMSD, `-9.398 kcal/mol`), directly validating why blind `FoldedCore_Box` screening cannot guarantee native binding-site identification.

### 4.3 Technical Feasibility of Retaining `HEM` in `2D0T`/`4PK5` with Existing Installed Software

Using existing WSL packages (`meeko 0.8.0`, `vina 1.2.7`, `rdkit 2026.09.2`) without `obabel`, `autogrid4`, `autodock4`, or `prepare_receptor4.py`:

1. **Meeko `0.8.0` lacks a `HEM` residue template:** `residue_chem_templates.json` (813 entries) contains isolated `FE` ion only, not `HEM`. Running `mk_prepare_receptor.py` on chain A of `2D0T` (`HEM A:404`) or `4PK5` (`HEM A:501`) with `HEM` present fails (`rc = 1`, both default and `--allow_bad_res -x` delete-bad-residues) with:
   - `RuntimeWarning: Input residues {'A:404': 'HEM'} not in residue_templates`
   - `Failed to create rdkitmol from cif. Error: Explicit valence for atom # 13 O, 2, is greater than permitted -> template for HEM will be None.`
   - `AttributeError: 'NoneType' object has no attribute 'parent'`
   Removing `HEM` restores `rc = 0` (`2D0T_A_stripped_noHEM`, `4PK5_A_stripped_noHEM`).
2. **PDB `HETATM` bond-order degradation:** Standalone `Chem.MolFromPDBFile` on `2D0T` `HEM` (PDB lacks bond orders/formal charges) saturates the porphyrin macrocycle and protonates carboxylate oxygens into gem-diols (`REMARK SMILES CCC1C(C)C2CC3C(C)C(CCC(O)O)...`), assigning `Fe` partial charge `+0.491` with atom type `Fe`.
3. **Vina `1.2.7` metal-typing and AD4 limitations:**
   - Uppercase PDBQT type `FE` causes a fatal Vina parser error (`Atom type FE is not a valid AutoDock type`).
   - Lowercase `Fe` (`[#26]`) is recognized by `vina 1.2.7` (`sf_name="vina"`) only as a generic donor sphere (`vina_ri=1.2, vina_donor=True`) lacking heme coordination chemistry, polarization, or redox electrostatics.
   - `Vina(sf_name="ad4").set_receptor(...)` aborts (`rc = 1`: `ERROR: Only flexible residues allowed with the AD4 scoring function. No (rigid) receptor`) because `AD4` grid maps require `autogrid4` (and standard AD4 receptor prep via `prepare_receptor4.py`/`obabel`), none of which are installed.

**Conclusion:** True heme-retaining holo-IDO1 docking is unsupported out-of-the-box by the installed `Meeko 0.8.0` + `Vina 1.2.7` stack; future holo-IDO1 modeling requires custom CCD-templated `HEM` parameterization (e.g., `AD4`/`autogrid4` maps or Vina-Zn metal-aware workflows), alongside experimental testing of whether weak indole acids displace heme like `A1IEE`.

---

## 5. Recommended Manuscript Revisions, Figure 6 Caption & Limitation Text

Recommended drop-in edits (also saved standalone in `02_results/JTM_upgrade/docking/figure6_caption_and_limitations.md`):

1. **Methods (line 77) corrections needed for accuracy:**
   - Replace *"OpenBabel/Meeko PDBQT conversion"* with *"RDKit (v2026.09.2; ETKDGv3 with MMFF94/UFF optimization from neutral PubChem SMILES) and Meeko (v0.8.0) PDBQT conversion"* (since `obabel` was not used).
   - Explicitly note that `9FOH` (IDO1) and `9J9D` (TGFBR1) were mmCIF-only and triggered automatic AlphaFold fallback; state exact search-box distinction (`26x26x26 A` `CoCrystal_Pocket` for `DAPK3` and `5NI9` Chain A vs. `FoldedCore_Box` `35,297–216,000 A^3` for the remaining 9 targets); and disclose the `5NI9` (`HLA-DRA` chain A) / `1JDH` (`CTNNB1` chain A) / `6BSC` (`MUC1` N-terminal SEA chain A) chain-selection caveats—or restrict interpretation away from `ENO1`/`TCF7L2` docking scores.
2. **Figure 6 Caption & Results framing:** Reframe *IDO1* × indole-3-carboxylic acid (`-6.68 kcal/mol`) strictly as **apo (heme-free) cavity engagement** overlapping the apo-IDO1 inhibitor `A1IEE` (`9FOH`, `0.32 A`) and incompatible with holo-heme occupancy (`0.95–1.01 A` clash with `2D0T`/`4PK5` heme), displaying two nearly isoenergetic binding orientations (`4.34–4.36 A` RMSD flip, `0.01–0.04 kcal/mol` gap). Frame *PEMT* × urolithin A (`-7.84 kcal/mol`) as an **unrefined in vacuo AlphaFold prediction** (contact pLDDT `87.94–96.25`, proximity `4.10 A` to catalytic `GLU180`, with a nearly isoenergetic `12.31 A`-shifted transmembrane decoy at `-7.76 kcal/mol`) lacking membrane lipid and SAM cofactor context.
3. **Discussion Limitations addition:** Add a dedicated structural-docking limitation item summarizing (i) AlphaFold fallback & heme omission in IDO1 (`0.95 A` clash with `2D0T` heme; Meeko `0.8.0` lacks `HEM` templates and `autogrid4` is unavailable), (ii) blind `FoldedCore_Box` search-volume inflation (`>27,000 A^3`, shown to fail crystal redocking of `9FOH` + `A1IEE` at `9.68 A` RMSD vs. `0.63–0.65 A` in focused `26 A` boxes), (iii) neutral protonation and rigid in vacuo receptor treatment lacking membrane lipids/SAM for PEMT, and (iv) PDB chain-selection artifacts (`5NI9` HLA-DRA, `1JDH` beta-catenin, `6BSC` split SEA monomer).
