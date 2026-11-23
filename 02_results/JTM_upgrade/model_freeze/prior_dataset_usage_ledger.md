# Prior Dataset Usage Ledger, Preprocessing Rules & Stage A Model Freeze Audit

**Audit Date (UTC):** `2026-11-22`  
**Frozen Model:** `GMSig-4` (`IFNG`, `IDO1`, `CXCL10`, `ENO1`)  
**Authoritative Model Object:** `02_results/GMSig4_Model_Objects.rds`  
**Canonical Scoring Function:** `00_scripts/26a_gmsig4_score.R` (`score_gmsig4`)  
**Freeze Script:** `00_scripts/26a_freeze_gmsig4.R` (12/12 checks passed; `verification_checks.csv`)

---

## 1. Frozen GMSig-4 Specification & Verification Summary

- **Training Dataset:** TCGA-BRCA Primary TNBC (`n = 172`, `OS_events = 29`, median follow-up `2.501` years, `time_unit = OS.time / 365.25`).
- **Model Fit:** Penalized Cox Ridge regression (`glmnet::cv.glmnet(z, y, family = "cox", alpha = 0, nfolds = 5)`, `seed = 2026`, `s = "lambda.min"`, `cox.ties = "breslow"` in `glmnet 5.0`, `glmnet_internal_standardize = TRUE`).
  - Reconstructed `lambda.min = 0.111047382457078` (`lambda.1se = 90.0739555894297`).
  - Maximum absolute coefficient reconstruction error vs `GMSig4_Model_Objects.rds$weights`: `0`.
- **Frozen Weights (`GMSig4_Model_Objects.rds$weights`):**
  - `IFNG`: `-0.251802755687135`
  - `IDO1`: `-0.150867350934661`
  - `CXCL10`: `-0.0719280374738463`
  - `ENO1`: `+0.184881748437269`
- **Training Cohort Gene-Level Statistics (`training_scaling.csv`):**
  - `IFNG`: mean `2.87903488372093`, sample SD (`n - 1`) `2.01750909927975`
  - `IDO1`: mean `8.32608895348837`, sample SD (`n - 1`) `2.50993038468383`
  - `CXCL10`: mean `10.0687930232558`, sample SD (`n - 1`) `2.14571055722566`
  - `ENO1`: mean `14.8409191860465`, sample SD (`n - 1`) `0.782323262505113`
  - Training raw score SD: `0.446545188201463`; training median: `-0.0102359731895073`.
- **All 12 Assertions in `verification_checks.csv` Passed (`TRUE`):**
  - TCGA score max abs error: `3.5527e-15` (`< 1e-12`)
  - TCGA z-expression max abs error: `5.1070e-15` (`< 1e-12`)
  - TCGA risk groups identical: `TRUE`
  - Ridge reconstruction & explicit fold allocation coefficient max abs error: `0` (`< 1e-12`)
  - Gene/weight ordering invariant: `0` (`< 1e-12`)
  - Fail-fast rejection checks (missing gene, non-finite input, zero gene SD, single sample): all `TRUE`
  - `GSE58812` & `GSE96058_TNBC_Ext` saved z-expression weighted sum max abs error: `5.1070e-15` and `5.3291e-15` (`< 1e-12`).

---

## 2. Exact Expression Preprocessing, Log-Scale & Multi-Probe Aggregation Rules

Extracted directly from `00_scripts/19a_prepare_cohort_matrices.R`, `00_scripts/19a2_prepare_gse96058.R`, `00_scripts/07`, `00_scripts/09`, `00_scripts/11`, `00_scripts/19b..19g`, `00_scripts/20a..20b`, and `00_scripts/26a_gmsig4_score.R`:

| Cohort | Source File(s) | Observed Scale & Log-Transformation Rule | Multi-Probe / Gene-Level Aggregation Rule | Within-Cohort Scaling & Scoring Rule |
|---|---|---|---|---|
| **TCGA-BRCA** (`19a`, `07`) | `TCGA_BRCA_HiSeqV2.tsv.gz` | Already `log2(RSEM + 1)` from UCSC Xena (`20,530 x 1,218`, range `[0.0000, 20.9784]`, 99th pct `13.3864`). **No additional log transform** applied. | None (already unique gene symbols). | Subsets primary TNBC (`ER- & PR- & HER2-`, `n = 172`, `OS.time > 0`). Gene-wise `scale()` (center=mean, scale=sample SD `n-1`), then `Z %*% w`. |
| **GSE58812** (`19a`, `11`) | `GSE58812_series_matrix.txt.gz` (GPL570 MAS 5.0) | Raw MAS 5.0 linear intensities (`54,675 x 107`, range `[0.0523, 74650.2600]`, 99th pct `12279.6860`). Rule: `if (max(raw58, na.rm = TRUE) > 50) raw58 <- log2(pmax(raw58, 0) + 1)` -> triggered (`max > 50`), yielding `[0.0898, 16.1879]`, 99th pct `13.2122`. | Mapped via `AnnotationDbi::mapIds(hgu133plus2.db, keys=rownames, column="SYMBOL", keytype="PROBEID", multiVals="first")`; dropped `NA`/empty; aggregated across probes per symbol by **arithmetic mean on log2 scale** -> `21,358` genes. | Gene-wise `scale()` across `n = 107` TNBC samples, then `Z %*% w`. |
| **GSE96058** (`19a2`) | `GSE96058_gene_expression_n3273.csv.gz` + `GPL11154`/`GPL18573` series matrices | Already author-transformed `log2(FPKM)` (`30,865 x 385` cached TNBC/Basal, range `[-3.3219, 23.0488]`, 99th pct `7.9262`). **No additional log transform** applied. | None (author-provided gene-level matrix). Replicates (`grepl("repl", sample_title, ignore.case=TRUE)`) excluded. | Independently `scale()`d within each evaluated subset (`Strict n=143`, `Ext n=242`, `Basal n=339`), then `Z %*% w`. |
| **GSE194040** (`19a`, `09`) | `GSE194040_geneLevel_n988.txt.gz` + `GPL30493`/`GPL20078` series matrices | Already log2-normalized (75th-quantile aligned + `9.5` offset) and ComBat batch-corrected by authors (`19,134 x 988`, range `[2.3354, 21.2500]`, 99th pct `14.4116`). **No additional log transform** applied. | None (authors mean-collapsed probes to `19,134` common genes). Sample IDs matched via `sub("^ISPY2_", "", title)` and `sub("^X", "", colnames)`. | **Subgroup-wise** `scale()` inside each evaluated arm/subtype (`Pembro_TNBC`, `Pembro_All`, `Pembro_HRp`, `Paclitaxel_Ctrl_TNBC`, `Paclitaxel_Ctrl_All`) in `19b`–`19g`, then `Z %*% w`. |
| **GSE173839** (`19a`, `11`) | `GSE173839_durva_probeLevel_n105.txt.gz` + `GSE194040_ProbeAnnotation.txt.gz` + `GSE173839_biomarkers.csv.gz` | Already Agendia `gNormalizedSignal` (`log2`) (`range [3.3910, 19.3970]`, 99th pct `13.7470`). **No additional log transform** applied. | Exact `match(durva$ProbeName, annot$ProbeName)` -> `annot$GeneName`; dropped `NA`, `""`, `"---"`; aggregated probes per symbol by **arithmetic mean on log2 scale** -> `16,947` genes. | **Subgroup-wise** `scale()` inside each evaluated arm/subtype (`DurvaOlap_All`, `DurvaOlap_TNBC`, `Ctrl_All`, `Ctrl_TNBC`) excluding `pCR.status == -1`, then `Z %*% w`. |
| **GSE176078** (`20a`, `20b`) | `GSE176078_Wu_etal_2021_BRCA_scRNASeq.tar.gz` -> `GSE176078_TNBC_annotated.h5ad` | Raw UMI counts subset to `subtype == "TNBC"` (`42,512` cells, `10` patients), normalized via `sc.pp.normalize_total(target_sum=1e4)` + `sc.pp.log1p` (`ln(CP10K + 1)`). | `adata.var_names_make_unique()`; patient pseudobulk = arithmetic mean of `ln(CP10K + 1)` across cells per patient (`orig.ident`). | Z-scored across the `10` TNBC patients (`ddof=1`); `20b` computed `Risk_C16_IfngIdo1Cxcl10Eno1` using **UniBeta** weights, not Ridge weights. |
| **GSE103091** (`06`, `08`) | `GSE103091_series_matrix.txt.gz` (GPL570 MAS 5.0) | Raw MAS 5.0 linear scale (`54,675 x 238`, range `[0.0523, 139282.6000]`, 99th pct `12315.0361`). **Never preprocessed, log-transformed, or cached.** | Never aggregated in prior pipeline. | Never scored in prior pipeline. |

### Scoring Invariant (`26a_gmsig4_score.R`)
`score_gmsig4(expr, weights, samples)` assumes `expr` is already on the cohort's documented log expression scale, applies **no automatic log transformation**, enforces finite values, unique gene row names, unique sample column names, `length(samples) >= 2`, and `sd > 0` for all 4 genes, performs gene-wise cohort centering and sample SD (`n - 1`) scaling via `scale(x, center = TRUE, scale = TRUE)`, computes `score = as.numeric(z %*% weights[genes])` (no post-sum score rescaling), and labels `RiskGroup = ifelse(score > median(score), "High", "Low")`.

---

## 3. Complete Prior Dataset Usage Ledger & Leakage Assessment

Full machine-readable table: [`prior_dataset_usage_ledger.csv`](file:///D:/Bioinformatics/TNBC_Microbiome/02_results/JTM_upgrade/model_freeze/prior_dataset_usage_ledger.csv).

### Summary of Transcriptomic & Single-Cell Cohorts

| Dataset | Primary Role | All Pipeline Roles | Clean Independent External Validation for GMSig-4? | Rationale / Evidence in `00_scripts/07..25` |
|---|---|---|---|---|
| **TCGA-BRCA** | `trained` | `trained; used_for_gene_selection; used_for_model_selection; evaluated` | **No** | Primary training set (`n=172` TNBC, 29 OS events) for DE/Cox gene screen (`07`, `19b`, `19c`), model selection (`19b`, `19d`, `19e`), Ridge Cox fit (`19f`), and nomogram (`19f`, `19h`, `25b`). |
| **GSE58812** | `used_for_model_selection` | `used_for_gene_selection; used_for_model_selection; evaluated` | **No** | Used in multi-cohort gene evidence ranking (`11`, `17`, `18`) and model ranking/selection (`19b` sorted by `-g58_OS_C`; `19d`, `19e` required `g58_OS_KM_p < 0.05`). |
| **GSE96058** | `used_for_model_selection` | `used_for_gene_selection; used_for_model_selection; evaluated` | **No** | `19c` screened all gutMGene genes in `Strict`/`Ext`/`Basal` (`TCGA_GSE58812_GSE96058_Univariate_Cox_AllGenes.csv`), where `ENO1` was strongly prognostic (`Strict p=2.62e-4`, `Ext p=1.49e-3`, `Basal p=1.82e-4`) despite being non-significant in TCGA (`p=0.272`) and GSE58812 (`p=0.276`). `19d`/`19e` included GSE96058 in 5-cohort selection (`Signature_Candidate_Selection_All5Cohorts.csv`). |
| **GSE194040** | `used_for_model_selection` | `used_for_gene_selection; used_for_model_selection; evaluated` | **No** | `09` screened gutMGene targets for Pembro pCR; `11`/`17`/`18` fed Pembro pCR into `evidence_score`; `19b`/`19d`/`19e` ranked candidate models on `pembro_TN_AUC`. |
| **GSE173839** | `used_for_model_selection` | `used_for_gene_selection; used_for_model_selection; evaluated` | **No** | `11`/`17`/`18` fed Durva pCR into `evidence_score`; `19b`–`19e` filtered/ranked candidates requiring presence in GSE173839 (`all5_ici == TRUE`, favoring `ENO1` over missing `TGFBR1`/`DAPK3`) and `durva_All_AUC`. |
| **GSE103091** | `downloaded_unused` | `downloaded_unused` | **No** (unused; overlaps GSE58812) | Downloaded (`06`, `06c`, `06d`) and metadata-inspected (`08`), but **never processed or evaluated** (see Section 4). |
| **GSE176078** | `evaluated` | `evaluated` | **No** (scRNA-seq mechanistic only) | Not used in gene/model selection or survival/ICI benchmarking; used downstream in `20a`/`20b`/`22_fig5`/`25c` for cell-type localization and pseudobulk TME/CellComm correlation. |
| **GSE169246** | `downloaded_unused` | `downloaded_unused` | **No** | Only an orphan barcode file (`GSE169246_TNBC_RNA.barcode.tsv.gz`) exists in `01_data/scRNA/`; zero references in `00..25`. |

**Critical Conclusion:** Because all five bulk cohorts (`TCGA-BRCA`, `GSE58812`, `GSE96058`, `GSE194040`, `GSE173839`) participated in gene screening and/or candidate-panel selection prior to locking GMSig-4 in `19f`, **none of the previously evaluated cohorts constitute an independent, leakage-free external validation of GMSig-4**.

### Summary of Microbiome, Knowledgebase, eQTL & GWAS Public Datasets
- **gutMGene v2.0 (`gutMGene_v2_human`)**: `used_for_gene_selection` (`00`, `01`, `02`, `05`, `07`, `19b`–`19d`, `20a`–`20b`, `21b`) — defined the 238-gene master universe.
- **GMrepo (`PRJNA658160` n=940, `PRJNA383849` n=288, `PRJNA872152` n=43)**: `used_for_gene_selection` (`00`, `01`, `03`, `04`, `05`, `22_fig1`) — differential BC genera in `04` anchored the 113-gene `GMrepo_anchored_gutMGene_targets.csv` in `05`.
- **HGMT 8 BC Cohorts (`PRJNA345373` n=94, `PRJNA383849` n=288, `PRJNA453965` n=133 WGS+KO, `PRJNA494824` n=31/8 BC WGS+KO, `PRJNA566060` n=36 BC, `PRJNA658160` n=935, `PRJNA718520` n=76 BC WGS+KO, `PRJNA751260` n=36)**: `evaluated` (`13`, `15`, `15b`, `16`, `22_fig1`, `23_audit`). Note BioProjects `PRJNA658160` and `PRJNA383849` overlap with GMrepo.
- **eQTLGen Blood cis-eQTL (`N = 31,684`)** & **BCAC `ieu-a-1126` (Overall), `ieu-a-1128` (ER-), `ieu-a-1135` (OncoArray ER-)**: `used_for_gene_selection; evaluated` (`12`, `14`, `14b`, `14f`, `14g`, `17`, `18`, `22_fig2`, `25d1`, `25d3`) — eQTLGen->BCAC MR in `17` fed `+1.5` into `evidence_score` in `18` (`used_for_gene_selection` for causal target panels, though none of `IFNG`, `IDO1`, `CXCL10`, `ENO1` were eQTLGen MR hits).
- **MiBioGen (`211` taxa)**, **GTEx V8 Breast (`N = 396`)**, and **Zhang et al. 2020 BCAC (`GCST010098`, `GCST010100`, `GCST90454344`)**: `evaluated` (`06*`, `12`, `13`, `14*`, `22_fig2`, `23*`, `25*`) — forward/reverse MR, Steiger, and coloc validation; not used in GMSig-4 scoring/training.

---

## 4. Clarification: Script `11` vs Scripts `19*` on `GSE103091`

1. **Download & Inspection (`06`, `06c`, `06d`, `08`):** `GSE103091_series_matrix.txt.gz` (`54,675` Affymetrix GPL570 probes x `238` samples, linear MAS 5.0 intensities `[0.0523, 139282.6000]`) was downloaded in `06`/`06c`/`06d` and inspected for dimensions/characteristics in `08_inspect_geo_cohorts.R`.
2. **Script `11` (`11_validate_gse58812_gse103091_gse173839.R`):** Despite carrying `gse103091` in its filename, **script `11` contains zero code referencing, loading, or processing `GSE103091`**. It loads only `GSE58812_series_matrix.txt.gz` (`n = 107`) and `GSE173839_durva_probeLevel_n105.txt.gz` (`n = 105`), merging them with `TCGA_TNBC_gutMGene_screen.csv` (`07`) and `ISPY2_Pembrolizumab_gutMGene_pCR_screen.csv` (`09`) into `MultiCohort_TNBC_gutMGene_evidence_matrix.csv`.
3. **Scripts `19*` (`19a`, `19a2`, `19b`–`19h`):** None of the `19*` scripts read, cache, or evaluate `GSE103091` (`01_data/TCGA_GEO/clean_cache/` contains no `GSE103091` object).
4. **Verified Overlap with `GSE58812`:** Direct comparison of `GSE58812_series_matrix.txt.gz` (`107` GSMs) and `GSE103091_series_matrix.txt.gz` (`238` GSMs) confirms that **all 107 GSM accessions and sample titles (`GP_001`..`GP_107`, labeled `diagnosis: triple negative breast cancer`) of `GSE58812` are 100% duplicated inside `GSE103091`**, alongside `131` non-overlapping samples (`ihc: TNBC`). Any future analysis using `GSE103091` must explicitly deduplicate these 107 GSMs against `GSE58812`.

---

## 5. Confirmed Stage A Limitations

1. **`Continuous_HR_perSD` Reporting Bug in `00_scripts/19f_signature_clinical_nomogram.R` (Lines 49 & 66):**
   - Line 49 fits `cx <- coxph(y ~ risk_score)` on the **unscaled** linear predictor `risk_score`, and line 66 outputs `cx_sum$conf.int[1, "exp(coef)"]` under the column name `Continuous_HR_perSD` in `02_results/GMSig4_Prognostic_Validation_Summary.csv`.
   - Because cohort `risk_score` SDs range from `0.4465` to `0.4836` (not `1.0`), `19f` reported the HR per **1.0 raw GMSig-4 unit** rather than per 1 SD (`coxph(y ~ scale(risk_score))`):
     - `TCGA_TNBC_OS` (`n=172`, `29` events, `score_sd=0.446545`): `19f` reported `4.3511` vs true per-SD `HR = 1.9282` (`95% CI 1.3023–2.8551`)
     - `GSE58812_TNBC_OS` (`n=107`, `29` events, `score_sd=0.458346`): `19f` reported `4.7190` vs true per-SD `HR = 2.0364` (`95% CI 1.3678–3.0317`)
     - `GSE58812_TNBC_MFS` (`n=107`, `31` events, `score_sd=0.458346`): `19f` reported `4.4174` vs true per-SD `HR = 1.9757` (`95% CI 1.3458–2.9003`)
     - `GSE96058_TNBC_Strict_OS` (`n=143`, `26` events, `score_sd=0.462193`): `19f` reported `4.3601` vs true per-SD `HR = 1.9750` (`95% CI 1.2368–3.1540`)
     - `GSE96058_TNBC_Ext_OS` (`n=242`, `50` events, `score_sd=0.468185`): `19f` reported `3.8751` vs true per-SD `HR = 1.8855` (`95% CI 1.3668–2.6009`)
     - `GSE96058_Basal_OS` (`n=339`, `67` events, `score_sd=0.483581`): `19f` reported `3.5308` vs true per-SD `HR = 1.8405` (`95% CI 1.4035–2.4135`)
   - Similarly, `19f` multivariable Cox models (`mv_tcga`, `mv_g960`, lines 116 & 140) fit raw unscaled `GMSig4` (`HR = 3.3088` in TCGA; `2.3173` in GSE96058_TNBC_Ext). Discrimination (`C_index`), log-rank `P`, Wald `P`, `HighLow_HR`, and time-dependent AUCs are unaffected by linear score scaling.
2. **Uncalibrated Baseline Risk for Standalone Ridge Model:**
   - `GMSig4_Model_Objects.rds` did not save the `cv.glmnet` object, CV fold assignments (`folds_originally_saved = FALSE`), or any standalone Ridge Cox baseline cumulative hazard $H_0(t)$ (`original_ridge_baseline_saved = FALSE`).
   - `nomo_tcga` and `nomo_g960` are separate unpenalized `coxph` models including clinical covariates (`Age`, `Stage_III_IV`/`Tumor_mm`, `Node_Pos`); their exported baselines (`nomo_*_existing_baseline_uncentered.csv`) cannot substitute for the Ridge baseline (`absolute_risk_validated = FALSE`).
3. **Cohort-Relative Z-Scaling & Single-Patient Non-Portability:**
   - All scoring (`19b`–`19g` and `26a_gmsig4_score.R`) relies on within-cohort (and in `19g`, within-subgroup/arm) gene-wise `scale()` (`mean = 0, sample SD = 1`), and `High`/`Low` groups use within-cohort `score > median(score)`.
   - Thus scores represent within-cohort relative rankings, single-sample inputs are rejected, and TCGA training means/SDs (`training_scaling.csv`) are not valid across external microarray/RNA-seq platforms without prespecified calibration (`single_sample_validated = FALSE`).
4. **Prior Bulk-Cohort Selection Circularity:**
   - Every bulk cohort previously evaluated (`TCGA`, `GSE58812`, `GSE96058`, `GSE194040`, `GSE173839`) participated in gene screening and/or model selection before locking `GMSig-4`.
5. **Downstream scRNA Weight Mismatch (`20b` / `25c`):**
   - `20b_scrna_tnbc_summary_and_tme.py` (and `scRNA_TNBC_PatientLevel_CellComm.csv` plotted in `25c`) scored the 10 `GSE176078` TNBC pseudobulks using C16 `UniBeta` weights (`-0.5941, -0.4378, -0.3292, +0.2195`) rather than the locked GMSig-4 Ridge weights (`-0.2518, -0.1509, -0.0719, +0.1849`).
6. **Historical Input Guards & Cox Tie Handling:**
   - Historical `scale()` lacked fail-fast checks for missing genes, duplicate gene symbols, non-finite expression, zero SD, or `n < 2` (now enforced by `score_gmsig4`).
   - `19f` used `glmnet 5.0` default `cox.ties = "breslow"` (documented in `model_manifest.json`; `glmnet 5.1` plans to switch default to `"efron"`).
