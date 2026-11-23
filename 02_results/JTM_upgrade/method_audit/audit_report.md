# Directed Methodological Audit Report (JTM Optimization Plan — Step 0)

**Audit Scope:** Read-only forensic audit of `D:\Bioinformatics\TNBC_Microbiome` and `D:\Bioinformatics\TNBC_Microbiome_GitHub` against `JTM_optimization_plan.md`, `04_manuscript\Manuscript_JTM.md`, and `04_manuscript\Additional_file_1_JTM.md`.  
**Deliverables Created (exclusively under `02_results/JTM_upgrade/method_audit/`):**
- `audit_report.md` (this report)
- `scrna_bh_recomputation.csv` (911 major lineage interactions with side-by-side buggy vs true BH `q_BH` verified across `scipy`, `R`, and manual step-up implementation)
- `manuscript_overclaims_table.csv` (22 structured audit entries with exact line numbers and exact replacement wording)
- `provenance_25c_scrna_cell_comm.py` (read-only exact copy of the generating script `D:\Bioinformatics\TNBC_Microbiome_GitHub\00_scripts\25c_scrna_cell_comm.py`)
- `source_provenance.json` (SHA256 hashes, byte sizes, timestamps, and git commit provenance)

---

## Executive Summary of Critical & High-Severity Findings

1. **Task 1 (`GSE173839` `STAT1_sig`, `PD1`, `PDL1`):**
   - `STAT1_sig` in `GSE173839_biomarkers.csv.gz` is **NOT** single-gene `STAT1`; it is a trial-provided multi-gene signature (`!Series_summary`: *"13 expression based signatures"*; Pearson `r = 0.8849` vs `STAT1` probe `A_24_P274270`, lower than its correlations with `CXCL10` `r = 0.9162`, `TAP1` `r = 0.9158`, and `CXCL11` `r = 0.9107`).
   - Meanwhile, `PD1` matches single-gene `PDCD1` (`A_23_P136405`, `r = 0.99999997`, `rho = 1.000000`), and `PDL1` correlates `r = 0.9628` with mean-collapsed `CD274` (`A_23_P256487`, `A_23_P338479`).
   - Although `19a_prepare_cohort_matrices.R` (L88–107) constructed gene-level expression `ex173` containing single-gene `STAT1`, `CD274`, and `PDCD1`, `19g_tme_immune_ici_analysis.R` (L270–272) and `24_preflight_checks.R` (L61–67) instead used `g173$clin$STAT1_sig`, `PD1`, and `PDL1` and labeled them `STAT1_AUC`, `PDCD1_AUC`, `CD274_AUC` (`GMSig4_Immunotherapy_pCR_Validation_Summary.csv`) and `STAT1`, `PDCD1`, `CD274` (`GMSig4_ICI_Subgroup_ExactScores.csv`).
   - Downstream, `Manuscript_JTM.md` explicitly equates `STAT1_sig` to single-gene `STAT1` (L69) and literally calls it `"single-gene ... STAT1_sig"` (L132), while `Additional_file_1_JTM.md` Table S5 (L114, L121–124) reports `STAT1_sig` AUCs (`0.694`, `0.758`, `0.682`, `0.747`) directly under `*STAT1* AUC` with no footnote.

2. **Task 2 (scRNA `25c_scrna_cell_comm.py` Provenance, Null Hypothesis & BH Bug):**
   - **Generating source located:** `D:\Bioinformatics\TNBC_Microbiome_GitHub\00_scripts\25c_scrna_cell_comm.py` (git commit `cf6975b7cdb21c487bc862e978c66222c2b81ba6`, SHA256 `2C7B19456FA98A2EAEC6A98E0631A0D902B7FF4C1D5906C788C7D9D2428B231F`), which generated the byte-identical scRNA CSVs in both `D:\Bioinformatics\TNBC_Microbiome\02_results\` and `D:\Bioinformatics\TNBC_Microbiome_GitHub\02_results\` (main `scRNA_TNBC_CellComm_Interactions.csv` SHA256 `F1C82B0F576793816B9D8CC8136B3B4D1D93A52936D540B0E25A48B5199C3490`). The copy at `D:\Bioinformatics\TNBC_Microbiome\00_scripts\25c_scrna_cell_comm.py` (SHA256 `B9F24441F832AE0806C5408D4B7964BE86056D3235C46E5347222B8FF6975C31`) was overwritten ~27 min after the GitHub commit into a Figure S3 plotting-only script.
   - **Wrong null hypothesis attribution:** `compute_comm()` permutes cell-type labels *within a single pool* (`1,000` perms overall across all 42,512 cells; `500` perms within Low and High pools; `300` perms within 12 sub-niches). Thus `perm_p` and `perm_fdr` test **overall pooled cell-type co-expression enrichment** against random cell-type labels—**NOT High-vs-Low GMSig-4 difference**. No pooled High-vs-Low permutation difference test exists in `25c`.
   - **Catastrophic Benjamini–Hochberg bug (`25c` L147–148):** `rk = out["perm_p"].rank(method="min"); out["perm_fdr"] = np.minimum(1.0, out["perm_p"] * len(out) / rk)` without step-up reverse cumulative minimum assigned `rk = 1` to all **321** ties at the minimum p-value (`perm_p = 1/1001 = 0.000999`), inflating their `perm_fdr` to `0.000999 * 911 / 1 = 0.910090` while assigning larger p-values (`perm_p = 0.001998`, `rk = 322`) `perm_fdr = 0.005653`.
   - **True BH recomputation (`scipy` / `R` / manual step-up, max difference `<= 2.81e-16`):** **331 / 911** interactions reach `q_BH < 0.05` (not `10 / 911`), including **`Myeloid CXCL10->CXCR3`** (`perm_p = 0.000999`, true `q_BH = 0.002835` vs buggy `0.910090`) and **`B-cells CXCL10->CXCR3`** (`perm_p = 0.004995`, true `q_BH = 0.013831` vs buggy `0.013916`). Claims in `Manuscript_JTM.md` (L137, L174, L374) and `Additional_file_1_JTM.md` (L50, L161–174) calling `B-cells CXCL10->CXCR3` a `"Significant pooled loss in High GMSig-4"` and rationalizing `Myeloid CXCL10->CXCR3 perm_fdr = 0.910` as `"due to active expression in both pools"` are completely invalidated.

3. **Task 3 (`14i` GTEx Coloc/MR & `19f` `Continuous_HR_perSD` Audit):**
   - **`19f` mislabeled continuous HR:** `19f_signature_clinical_nomogram.R` (L49, L66) fits `coxph(y ~ risk_score)` **without `scale(risk_score)`** (`SD = 0.4465–0.4836`), reporting raw-unit HRs (`3.53–4.72`) under column `Continuous_HR_perSD` in `GMSig4_Prognostic_Validation_Summary.csv`, whereas true per-1-SD HRs are **`1.841–2.036`** (while `SingleGene_Checkpoint_Survival_Benchmark.csv` in `19g` L324 *did* standardize to SD=1). Similarly, `19f` multivariable Cox models (L116, L140) report raw-unit HRs (`3.31` TCGA, `2.32` GSE96058) without stating the unit (`per-SD` equivalents: **`1.72 [1.15–2.57]`** and **`1.49 [1.08–2.05]`**).
   - **`14i` GTEx coloc gating & MR FDR:**
     - `PP.H3 + PP.H4 >= 0.50` is presented as a `"validity gate"` / `"required"` threshold (`Manuscript_JTM.md` L60, L100; `Additional_file_1_JTM.md` L25) rather than a study-chosen heuristic filter, omits `coloc.abf`'s single-causal-variant and prior-sensitivity assumptions, and was omitted in `17_eQTLGen_OpenGWAS_SMR_coloc.R` (applied only during plotting in `22_fig2`).
     - `Manuscript_JTM.md` L156 falsely claims **only** `LDHA` and `ST6GAL1` reached `PP.H3 + PP.H4 >= 0.50` in eQTLGen, whereas **9 / 17** eQTLGen rows pass `>= 0.50` (7 driven by `PP.H3 >= 0.928`).
     - `Manuscript_JTM.md` L174 mislabels `PP.H0 = 0.943` at `PEMT` as a `"null prior"` instead of a **posterior probability**.
     - `14i` L166 partitions GTEx Breast MR BH correction **`by = method`** (13 IVW vs 6 Wald ratio), yielding `AKT1 FDR = 0.0831` instead of **`0.1214`** across all 19 tested genes (and producing `"13–19 tested"` in `Manuscript_JTM.md` L58).

---

## Task 1: Exact Trace of `STAT1_sig`, `PD1`, and `PDL1` in `GSE173839`

### 1.1 Raw Source Definitions (`01_data/TCGA_GEO/GSE173839*`)
- **`GSE173839_series_matrix.txt.gz`:**
  - `!Series_title`: `"Durvalumab with olaparib and paclitaxel for high-risk HER2-negative stage II/III breast cancer: Results from the adaptively randomized I-SPY2 platform trial"` (PMID `34143979`).
  - `!Series_summary` (L7): *"In this analysis, we assessed **13 expression based signatures** related to immune, DRD, proliferation and estrogen signaling as specific predictor of pathologic complete response (pCR) to durvalumab/olaparib."*
  - `!Series_overall_design` (L8): GPL20078 Agilent 44K arrays on FFPE biopsies (`N = 71` durvalumab/olaparib; `N = 34` control). Five control samples have `pcr: -1` (2 HR=0, 3 HR=1), leaving `n = 71` treated (21 TNBC/HR=0, 50 HR=1) and `n = 29` evaluable controls (17 TNBC/HR=0, 12 HR=1) with `pcr` in `{0, 1}`.
- **`GSE173839_biomarkers.csv.gz` (`105 × 20`):**
  - Contains 14 pre-computed standardized biomarker columns (all near mean `-0.010` to `-0.013`, SD `0.995` to `1.000` across `n = 105`): `PD1`, `PDL1`, `T.cells_sig`, `B.cells_sig`, `Dendritic.cells_sig`, `Mast.cells_sig`, `CD68`, `TIS_sig`, `STAT1_sig`, `TAMsurr_TcClassII_ratio_sig`, `PARPi7_sig.`, `Mitotic_sig`, `ESR1_PGR_ave`, `SET.index`.

### 1.2 Empirical Verification Against `GSE173839_durva_probeLevel_n105.txt.gz` & `19a` Gene Matrix
Using `19a_prepare_cohort_matrices.R` (L88–107) probe-to-gene mapping (`GSE194040_ProbeAnnotation.txt.gz` on `GSE173839_durva_probeLevel_n105.txt.gz`, verified identical at shared probes to `GSE173839_series_matrix.txt.gz`):

| Biomarker Column (`GSE173839_biomarkers.csv.gz`) | Mapped Gene & Probes in `durva_probeLevel_n105` | Pearson `r` vs Mean Gene Expr (`n=105`) | Spearman `rho` (`n=105`) | Max `|Z(gene) - bio|` | Identity Verification |
|---|---|---:|---:|---:|---|
| `PD1` | `PDCD1` (1 probe: `A_23_P136405`) | `0.99999997` (`1.000000`) | `1.000000` | `0.0114` (lin. resid `0.0008`) | **Yes:** Standardized single-gene `PDCD1` expression (`≈ (PDCD1 - 7.4217)/1.0505`). |
| `PDL1` | `CD274` (2 probes: `A_23_P256487`, `A_23_P338479`) | `0.962779` (probes `0.7991`, `0.8449`) | `0.960483` | `0.7524` (OLS `R²=0.9272`) | **Approximate:** Closely tracks `CD274` probes, not identical to simple unweighted 2-probe mean. |
| `STAT1_sig` | `STAT1` (1 probe: `A_24_P274270`) | `0.884933` | `0.870207` | `1.5067` | **No — Multi-gene Signature:** Correlates *more strongly* with `CXCL10` (`r=0.91620`), `TAP1` (`r=0.91575`), and `CXCL11` (`r=0.91073`) than with `STAT1` (`r=0.88493`). |

### 1.3 Downstream Mislabeling Trace & Exact AUC Divergence (`19a` -> `19g` / `24` -> CSVs / Figures / Manuscripts)
`19a_prepare_cohort_matrices.R` saved `list(expr = ex173, clin = bio)` to `GSE173839_clean.rds`. Thus `g173$expr` *does* contain single-gene `STAT1`, `CD274`, and `PDCD1`. However, `19g_tme_immune_ici_analysis.R` (L270–272) and `24_preflight_checks.R` (L61–67) evaluated `g173$clin$PDL1`, `g173$clin$PD1`, and `g173$clin$STAT1_sig` and exported them under `CD274_AUC`, `PDCD1_AUC`, and `STAT1_AUC` (`GMSig4_Immunotherapy_pCR_Validation_Summary.csv`) and `CD274`, `PDCD1`, `STAT1` (`GMSig4_ICI_Subgroup_ExactScores.csv`).

| GSE173839 Subgroup (`n`, `pCR`) | `STAT1_sig` AUC (Reported as `STAT1_AUC`) | True Single-Gene `STAT1` AUC (`ex173`) | Delta (`sig - gene`) | `PDL1` AUC (Reported as `CD274_AUC`) | True Single-Gene `CD274` AUC (`ex173`) | `PD1` & `PDCD1` AUC |
|---|---:|---:|---:|---:|---:|---:|
| `DurvaOlap_All` (`71`, `29`) | **0.757800** | 0.696223 | **+0.061576** | 0.649425 | 0.669951 | 0.682266 |
| `DurvaOlap_TNBC` (`21`, `9`) | **0.694444** | 0.712963 | **-0.018519** | 0.555556 | 0.592593 | 0.527778 |
| `Ctrl_All` (`29`, `7`) | **0.746753** | 0.675325 | **+0.071429** | 0.792208 | 0.779221 | 0.746753 |
| `Ctrl_TNBC` (`17`, `6`) | **0.681818** | 0.712121 | **-0.030303** | 0.712121 | 0.712121 | 0.712121 |

- **Impact on Figures:** `22_fig5_ici_and_scrna.py` L180 (`Figure5C`) plots `"CD274 (PD-L1) Durva (0.556)"` using `GMSig4_PatientScores_GSE173839_Durva.csv$PDL1` rather than single-gene `CD274` (`0.593`). `STAT1_sig` was not plotted in Figure 5C (which plots Pembro TNBC single-gene `STAT1`, `0.863`).
- **Impact on Manuscripts:** `Manuscript_JTM.md` L69 and L132 explicitly mislabel `STAT1_sig` as a `"single-gene"` comparator, and `Additional_file_1_JTM.md` Table S5 (L114, L121–124) lists `STAT1_sig` and `PDL1` AUCs under `*STAT1* AUC` and `*CD274* AUC` with zero distinction from GSE194040 single-gene AUCs.

---

## Task 2: scRNA Permutation Provenance, Null Hypothesis & True BH Recomputation

### 2.1 Source Provenance & File Integrity
- **Original Generating Script (copied read-only to `provenance_25c_scrna_cell_comm.py`):**
  - Path: `D:\Bioinformatics\TNBC_Microbiome_GitHub\00_scripts\25c_scrna_cell_comm.py` (Git commit `cf6975b7cdb21c487bc862e978c66222c2b81ba6`, `2026-11-20 21:52:21 +0800`)
  - Size: `20,900` bytes
  - SHA256: `2C7B19456FA98A2EAEC6A98E0631A0D902B7FF4C1D5906C788C7D9D2428B231F`
- **Overwritten Main Workspace Script (plotting-only script, NOT generating code):**
  - Path: `D:\Bioinformatics\TNBC_Microbiome\00_scripts\25c_scrna_cell_comm.py` (`2026-11-20 22:09:25`)
  - Size: `12,479` bytes
  - SHA256: `B9F24441F832AE0806C5408D4B7964BE86056D3235C46E5347222B8FF6975C31`
- **Generated CSV Outputs (byte-identical between main `02_results/` and GitHub `02_results/`):**
  - `scRNA_TNBC_CellComm_Interactions.csv` (`164,252` bytes, SHA256 `F1C82B0F576793816B9D8CC8136B3B4D1D93A52936D540B0E25A48B5199C3490`) — **Original raw `perm_p` column IS present for all `m = 911` pairs.**
  - `scRNA_TNBC_CellComm_HighVsLow_GMSig4.csv` (`164,358` bytes, SHA256 `C5AA8B633327230E64B8CFEC555B113E16C7BB5C916A15319EC151F58B30EDA4`)
  - `scRNA_TNBC_SubNiche_CellComm_HighVsLow.csv` (`277,800` bytes, SHA256 `7079E69576C6F683E870E5297B89975992D29A05D0C0802DF69C231C60C4448B`)

### 2.2 Exact Null Hypothesis & Permutation Counts
In `compute_comm(df_in, ct_col, lr_list, n_perm, seed)` (GitHub `25c` L105–149):
1. Active pairs tested satisfy `mean_lig > 0.02`, `pct_lig >= 0.03`, `mean_rec > 0.02`, `pct_rec >= 0.03` (where `pct = mean(expr > 0.1)`), with score `S = sqrt(mean_lig * mean_rec)`.
2. For each iteration `p = 1..n_perm`, `shuf = rng.permutation(n_cells)` permutes **cell-type labels across cells within `df_in`**, preserving cell-type sizes, and computes `null_scores = sqrt(null_mean_lig * null_mean_rec)` with `perm_p = (count(null_scores >= S) + 1) / (n_perm + 1)`.
3. **Three separate within-pool enrichment calls exist:**
   - **Overall major lineages (`comm_all`, L152):** `n_perm = 1000`, `seed = 42` across all `42,512` cells (`911` active pairs, min `perm_p = 1/1001 = 0.000999`). Only this call writes `perm_fdr` (`scRNA_TNBC_CellComm_Interactions.csv`) and merges `perm_fdr` onto `scRNA_TNBC_CellComm_HighVsLow_GMSig4.csv` (L159).
   - **High vs Low major pools (`comm_low`, `comm_high`, L157–158):** `n_perm = 500` each (`seed = 101` Low, `202` High; min `perm_p = 1/501 = 0.001996`), yielding `p_Low` and `p_High` (no BH adjustment, and **no test comparing High vs Low**).
   - **Sub-niche pools (`comm_sub_low`, `comm_sub_high`, L215–216):** `n_perm = 300` each (`seed = 303` Low, `404` High across 12 kept sub-niches; `2,102` outer-merged pairs, min `perm_p = 1/301 = 0.003322`), again testing within-pool enrichment only, with **no High-vs-Low difference test**.
- **Conclusion:** `perm_fdr` tests **pooled overall cell-type co-expression enrichment across all 10 tumors combined**, **NOT** High-vs-Low GMSig-4 difference. The only High-vs-Low hypothesis tests in the entire scRNA pipeline are the exploratory `n = 5 vs 5` patient-level Wilcoxon/Spearman tests (`scRNA_TNBC_PatientLevel_CellComm_Summary.csv`), all `p >= 0.117`.

### 2.3 Verification of the BH Math Bug & Exact Side-by-Side Recomputation
Lines 147–148 of GitHub `25c` execute:
```python
rk = out["perm_p"].rank(method="min")
out["perm_fdr"] = np.minimum(1.0, out["perm_p"] * len(out) / rk)
```
Because **321** pairs share the floor p-value `perm_p = 1 / 1001 = 0.000999001`, `rank(method="min")` assigns `rk = 1` to all 321, giving `0.000999001 * 911 / 1 = 0.910090` (reproduced in `scrna_bh_recomputation.csv` to `max|diff| = 9.03e-14`), while the 4 pairs with `perm_p = 0.001998` get `rk = 322` (`0.001998 * 911 / 322 = 0.005653`), inverting p-value monotonicity!

Recomputing true step-up Benjamini–Hochberg `q_BH` on exact `perm_p` (`m = 911`) via `scipy.stats.false_discovery_control(..., method='bh')`, `R stats::p.adjust(..., method='BH')`, and manual reverse cumulative minimum (`max|scipy - R| = 2.81e-16`, `max|scipy - manual| = 1.11e-16`):

| Metric / Threshold (`m = 911` Overall Major Lineage Pairs) | Buggy `25c` (`perm_fdr`) | True BH (`q_BH`) | Difference |
|---|---:|---:|---:|
| Nominal `perm_p < 0.05` | 336 | 336 | — |
| FDR `< 0.05` | **10** | **331** | **+321 pairs recovered** (0 lost) |
| FDR `< 0.01` | **4** | **325** | **+321 pairs recovered** |
| FDR `< 0.10` | **13** | **334** | **+321 pairs recovered** |
| Pairs with floor `perm_p = 0.000999` (`n = 321`) | `0.910090` (NS) | **`0.002835`** (Significant) | `-0.907255` |

**Key Biological & Manuscript-Cited Axes (Side-by-Side):**

| Axis (`source -> target | lr_pair`) | `comm_score` (Low / High) | `p_Low` / `p_High` (500-perm) | Overall `perm_p` (1000-perm) | Buggy `perm_fdr` | True `q_BH` | Status & Interpretation Fix |
|---|---:|---:|---:|---:|---:|---|
| `Myeloid -> T-cells | CXCL10->CXCR3` | 0.739 (0.738 / 0.564) | 0.0020 / 0.0020 | 0.000999 | **0.910090** | **0.002835** | Flipped NS -> Sig overall enrichment; `0.910` was pure bug |
| `Myeloid -> T-cells | CXCL9->CXCR3` | 0.608 (0.596 / 0.522) | 0.0020 / 0.0020 | 0.000999 | **0.910090** | **0.002835** | Flipped NS -> Sig overall enrichment |
| `Myeloid -> T-cells | CD274->PDCD1` | 0.154 (0.147 / 0.133) | 0.0020 / 0.0020 | 0.000999 | **0.910090** | **0.002835** | Flipped NS -> Sig overall enrichment |
| `T-cells -> Myeloid | IFNG->IFNGR1` | 0.584 (0.691 / 0.401) | 0.0020 / 0.0020 | 0.000999 | **0.910090** | **0.002835** | Flipped NS -> Sig overall enrichment |
| `T-cells -> Myeloid | IFNG->IFNGR2` | 0.790 (0.982 / 0.450) | 0.0020 / 0.0020 | 0.000999 | **0.910090** | **0.002835** | Flipped NS -> Sig overall enrichment |
| `CAFs -> T-cells | COL1A1->CD44` | 1.899 (1.931 / 1.860) | 0.0020 / 0.0020 | 0.000999 | **0.910090** | **0.002835** | Flipped NS -> Sig overall enrichment |
| `PVL -> B-cells | CXCL9->CXCR3` | 0.226 (0.248 / 0.138) | 0.0818 / 0.1557 | 0.001998 | 0.005653 | **0.005601** | Tie-rank minor shift; overall enrichment only |
| `T-cells -> T-cells | LGALS9->HAVCR2` | 0.185 (0.206 / 0.146) | 1.0000 / 0.0020 | 0.001998 | 0.005653 | **0.005601** | Tie-rank minor shift; overall enrichment only |
| `Myeloid -> PVL | TGFB1->TGFBR2` | 0.325 (0.366 / 0.236) | 0.6926 / 0.0439 | 0.001998 | 0.005653 | **0.005601** | Tie-rank minor shift; overall enrichment only |
| `CAFs -> Myeloid | ICAM1->ITGAL` | 0.197 (0.228 / 0.144) | 0.3872 / 0.0080 | 0.001998 | 0.005653 | **0.005601** | Tie-rank minor shift; overall enrichment only |
| `CAFs -> T-cells | LGALS9->HAVCR2` | 0.192 (0.236 / 0.167) | 0.4491 / 0.0020 | 0.003996 | 0.011167 | **0.011167** | Overall enrichment; NOT a High-vs-Low `"shift"` test |
| `B-cells -> T-cells | CXCL10->CXCR3` | 0.266 (0.286 / 0.000) | 0.5269 / 1.0000 | 0.004995 | 0.013916 | **0.013831** | Overall enrichment; NOT `"Significant pooled loss in High"` |
| `CAFs -> Myeloid | TGFB1->TGFBR2` | 0.316 (0.343 / 0.239) | 0.9960 / 0.0020 | 0.004995 | 0.013916 | **0.013831** | Overall enrichment; NOT a High-vs-Low `"shift"` test |
| `Endothelial -> CAFs | TGFB1->TGFBR2` | 0.327 (0.408 / 0.256) | 0.0180 / 0.0040 | 0.004995 | 0.013916 | **0.013831** | Tie-rank minor shift; overall enrichment only |
| `Cancer Epithelial -> Plasmablasts | CCL2->CCR2` | 0.103 (0.134 / 0.073) | 0.0160 / 0.0100 | 0.011988 | 0.033094 | **0.033094** | Overall enrichment only |
| `Endothelial -> Endothelial | MIF->CD74` | 2.063 (2.147 / 1.920) | 0.2275 / 0.0240 | 0.017982 | 0.049491 | **0.049491** | Overall enrichment only |

---

## Task 3: Audit of `14i` (GTEx Coloc/MR) & `19f` (`Continuous_HR_perSD` / Multivariable Scale)

### 3.1 `14i_gtex_breast_mr_and_tnbc_coloc.R` Audit
1. **Heuristic `PP.H3 + PP.H4 >= 0.50` Gating (`14i` L89–93):**
   - Technically prevents division of near-zero posteriors (`PEMT` `cond_H4_raw = 0.7039` when `PP.H0 = 0.9426`; `CLDN7` `0.6453` when `PP.H1 = 0.7214`; `VEGFA` `0.5090` when `PP.H0 = 0.8894`). Across 27 GTEx Breast loci, 26 fail `>= 0.50` (`cond_H4 = NA`) and only `LDHA` passes (`PP.H3 = 0.7390`, `PP.H4 = 0.00301`, `PP.H3 + PP.H4 = 0.7420`, `cond_H4 = 0.00406`).
   - **Audit findings:**
     - `Manuscript_JTM.md` L60, L100 and `Additional_file_1_JTM.md` L25 elevate this heuristic rule to `"substantial combined causal signal evidence"`, `"required for conditional colocalization interpretation"`, and a `"strict ... validity gate"` without stating that it is a study-chosen threshold or noting `coloc.abf`'s single-causal-variant assumption and prior sensitivity (`p1=1e-4, p2=1e-4, p12=1e-5`).
     - `17_eQTLGen_OpenGWAS_SMR_coloc.R` does not compute `PP.H3 + PP.H4` or `cond_H4` in `Coloc_eQTLGen_BCAC_results.csv`; gating was applied post-hoc in `22_fig2_mr_coloc_mediation.R` L315.
     - `Manuscript_JTM.md` L156 states falsely that **only** `LDHA` and `ST6GAL1` reached `PP.H3 + PP.H4 >= 0.50` in eQTLGen, whereas **9 / 17** eQTLGen tests reach `>= 0.50` (the other 7 having `PP.H3 = 0.928–1.000`).
     - `Manuscript_JTM.md` L174 mislabels `PP.H0 = 0.943` at `PEMT` as `"null prior"` rather than **null posterior probability**.
     - `Manuscript_JTM.md` L21 calls all mammary colocalization `"underpowered"`, overlooking 9 loci with strong single-trait or distinct-variant posterior mass (`LDHA` `PP.H3=0.739`; 6 `PP.H1 >= 0.691` eQTL-only loci including `AKT1` `0.786`; 2 `PP.H2 >= 0.549` GWAS-only loci `TCF7L2`, `ENO1`).
     - (`14i` L68 also hardcodes `s = 0.15` instead of `11488/118987 = 0.0965`, verified numerically invariant in `coloc.abf` when `beta` and `varbeta` are supplied.)
2. **GTEx Breast MR BH FDR Grouping (`14i` L166):**
   - `GTEx_Breast_clumped_IVs.csv` holds `20` genes (`F11R` drops in harmonization), yielding `19` tested genes in `MR_GTExBreast_GCST010100_TNBC_results.csv` (`13` IVW + `6` Wald ratio).
   - `14i` L166 runs `p.adjust(pval, method = "BH"), by = method`, splitting the 13 IVW genes from the 6 Wald-ratio genes. This reports `AKT1` `FDR = 0.08306` (`0.00639 * 13 / 1`) in `Manuscript_JTM.md` L98, `Additional_file_1_JTM.md` L23 & L94, and `Figure 2E`, instead of **`0.12140`** (`0.00639 * 19 / 1`) across all 19 tested GTEx Breast genes, and created the unclear `"13–19 tested"` wording in `Manuscript_JTM.md` L58.

### 3.2 `19f_signature_clinical_nomogram.R` Continuous & Multivariable HR Estimand Audit
`GMSig-4` (`IFNG -0.251803`, `IDO1 -0.150867`, `CXCL10 -0.071928`, `ENO1 +0.184882` on within-cohort gene Z-scores) has cohort score **SD = 0.4465 to 0.4836**, NOT `1.0`. Because `19f` fits Cox models on the unscaled score (`coxph(y ~ risk_score)` at L49, `mv_tcga` at L116, `mv_g960` at L140), every continuous GMSig-4 HR in `19f` is per **1 raw score unit (~2.1–2.2 SDs)**:

| Cohort / Model (`19f`) | Score SD | Reported Raw-Unit HR (95% CI) | Label in CSV / Manuscript | True Per-1-SD HR (95% CI) | Cox `p` |
|---|---:|---:|---|---:|---:|
| `TCGA_TNBC_OS` (univariable, `n=172`) | 0.4465 | **4.3511** (1.8066–10.4795) | `Continuous_HR_perSD` (`GMSig4_Prognostic...csv`) | **1.9282** (1.3023–2.8551) | 0.00104 |
| `GSE58812_TNBC_OS` (univariable, `n=107`) | 0.4583 | **4.7190** (1.9804–11.2443) | `Continuous_HR_perSD` (`GMSig4_Prognostic...csv`) | **2.0364** (1.3678–3.0317) | 0.00046 |
| `GSE58812_TNBC_MFS` (univariable, `n=107`) | 0.4583 | **4.4174** (1.9116–10.2082) | `Continuous_HR_perSD` (`GMSig4_Prognostic...csv`) | **1.9757** (1.3458–2.9003) | 0.00051 |
| `GSE96058_TNBC_Strict_OS` (univariable, `n=143`) | 0.4622 | **4.3601** (1.5837–12.0039) | `Continuous_HR_perSD` (`GMSig4_Prognostic...csv`) | **1.9750** (1.2368–3.1540) | 0.00438 |
| `GSE96058_TNBC_Ext_OS` (univariable, `n=242`) | 0.4682 | **3.8751** (1.9493–7.7033) | `Continuous_HR_perSD` (`GMSig4_Prognostic...csv`) | **1.8855** (1.3668–2.6009) | 0.00011 |
| `GSE96058_Basal_OS` (univariable, `n=339`) | 0.4836 | **3.5308** (2.0158–6.1843) | `Continuous_HR_perSD` (`GMSig4_Prognostic...csv`) | **1.8405** (1.4035–2.4135) | 1.03e-5 |
| `TCGA_TNBC` Multivariable (`n=168`) | 0.4511 | **3.3088** (1.3501–8.1090) | `"continuous GMSig-4"` (`Manuscript` L118, L356; `Table S4` L104; `Fig 3H`) | **1.7156** (1.1450–2.5707) | 0.00889 |
| `GSE96058_Ext` Multivariable (`n=234`) | 0.4714 | **2.3173** (1.1659–4.6059) | `"continuous GMSig-4"` (`Manuscript` L118, L356; `Table S4` L108; `Fig 3H`) | **1.4861** (1.0750–2.0543) | 0.01649 |

By contrast, `SingleGene_Checkpoint_Survival_Benchmark.csv` (`19g` L324) *did* call `scale()` before Cox fitting, so its `Continuous_HR_perSD` column is genuinely per 1 SD (with `STAT1` and `CD274` continuous Cox `p < 0.05` in **5/6** endpoints, whereas median-split KM log-rank was `< 0.05` in **2/6**).

### 3.3 Additional Overclaims Identified During Line-by-Line Audit
- **`Manuscript_JTM.md` L104:** Claims `"26 models retaining significance across all five core survival criteria and all ICI criteria"` and that `GMSig-4` `"ranked highest among Ridge candidates for pembrolizumab-TNBC discrimination"`. In `Signature_Candidate_Selection_All5Cohorts.csv` (`19e` L138–139), `26` is `all_surv_sig` (4 KM endpoints) AND `all5_ici` (`!is.na(durva_All_AUC)` = **array gene coverage in GSE173839**, NOT ICI significance; only `1/50` models had `p < 0.05` across all 4 ICI tests), and 3 Ridge models missing GSE173839 coverage (`C19` 0.847, `C2` 0.821, `C23` 0.784) had higher `pembro_TN_AUC` than `C16_Ridge` (`0.737`).
- **`Manuscript_JTM.md` L23, L184:** Uses `"immune-excluded phenotype"` / `"CAF-excluded tumor ecosystem"` without spatial transcriptomic or histological margin-vs-core data (and with non-significant patient-level CAF communication trends).

See `manuscript_overclaims_table.csv` for the complete 22-row line-by-line audit matrix with exact replacement text.
