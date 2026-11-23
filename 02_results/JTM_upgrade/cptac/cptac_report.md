# Stage C2: Narrow CPTAC Proteomics Orthogonal Validation (`brca_cptac_2020`)

**Date:** 2026-11-22  
**Cohort:** CPTAC-3 Prospective Breast Cancer (`brca_cptac_2020`, Krug et al., *Cell* 2020; PMID `33212010`; `N = 122` BRCA; **`n = 28` confirmed TNBC**, **`n = 16` strict clinical IHC TNBC**; **0 / 122 overlap** with TCGA-BRCA)  
**Frozen Model Authority:** [`frozen_coefficients.json`](file:///D:/Bioinformatics/TNBC_Microbiome/02_results/JTM_upgrade/model_freeze/frozen_coefficients.json) (`IFNG = -0.251802755687135`, `IDO1 = -0.150867350934661`, `CXCL10 = -0.0719280374738463`, `ENO1 = +0.184881748437269`)  
**Reproducible Scripts:** [`26h_download_cptac.py`](file:///D:/Bioinformatics/TNBC_Microbiome/00_scripts/26h_download_cptac.py), [`26h_cptac_proteomics_validation.py`](file:///D:/Bioinformatics/TNBC_Microbiome/00_scripts/26h_cptac_proteomics_validation.py)  
**Machine-Readable Outputs:** [`cptac_download_manifest.csv`](file:///D:/Bioinformatics/TNBC_Microbiome/02_results/JTM_upgrade/cptac/cptac_download_manifest.csv), [`cptac_tnbc_patient_table.csv`](file:///D:/Bioinformatics/TNBC_Microbiome/02_results/JTM_upgrade/cptac/cptac_tnbc_patient_table.csv), [`cptac_protein_concordance_summary.csv`](file:///D:/Bioinformatics/TNBC_Microbiome/02_results/JTM_upgrade/cptac/cptac_protein_concordance_summary.csv)

---

## 1. Executive Verdict & Scope Guardrails

1. **No 4-Protein GMSig-4 Score & Zero Imputation:** Global TMT-10-plex mass spectrometry (`data_protein_quantification.txt`, 9,715 proteins) detects **`IFNG` in 0 / 28 TNBC (0 / 122 BRCA)** and **`CD274` (PD-L1) in 0 / 28 TNBC (0 / 122 BRCA)**. Per protocol, no imputation was performed and no 4-protein GMSig-4 score was constructed. Survival endpoints (OS/DFS/RFS) are absent in `brca_cptac_2020`.
2. **Strong RNA-to-Protein Fidelity for Immune Targets, Modest for `ENO1`:** In non-imputed TNBC pairs, `CXCL10` ($\rho = 0.8805$, $p = 1.40 \times 10^{-7}$, $n = 21$), `IDO1` ($\rho = 0.8389$, $p = 5.70 \times 10^{-7}$, $n = 23$), and `STAT1` ($\rho = 0.8714$, $p = 1.60 \times 10^{-9}$, $n = 28$) show high RNA–protein concordance; `ENO1` shows modest positive concordance ($\rho = 0.4625$, $p = 0.0132$, $n = 28$). Missingness in `IDO1` protein (5 NA) is abundance-driven MNAR (median `RNA_IDO1` `0.185` in protein-NA vs `4.583` in detected, Mann–Whitney $p = 0.0128$), confirming that imputation is inappropriate.
3. **Orthogonal Protein-Level Validation of Frozen RNA `GMSig-4` Immune Axis:** Higher frozen RNA `GMSig-4` (poor-prognosis direction) strongly and inversely tracks non-imputed `IDO1` ($\rho = -0.7253$, $p = 9.01 \times 10^{-5}$), `CXCL10` ($\rho = -0.7740$, $p = 3.84 \times 10^{-5}$), and **disjoint orthogonal IFN/cytotoxic effector proteins** `STAT1` ($\rho = -0.8566$, $p = 6.03 \times 10^{-9}$), `GZMA` ($\rho = -0.7192$, $p = 1.61 \times 10^{-5}$), `CD8A` ($\rho = -0.6490$, $p = 2.50 \times 10^{-4}$), `PRF1` ($\rho = -0.6047$, $p = 0.00287$), and the 4-protein immune proxy ($\rho = -0.8626$, $p = 3.58 \times 10^{-9}$; High vs Low Wilcoxon $p = 1.50 \times 10^{-4}$), alongside `ESTIMATE_IMMUNE_SCORE` ($\rho = -0.7997$, $p = 3.30 \times 10^{-7}$; Wilcoxon $p = 4.78 \times 10^{-5}$).
4. **Honest Nulls — `CD3_TILS_COUNTS` (`n = 13`) and Composite `GMSig-4` vs Glycolysis Proteins (`n = 28`):** Pathologist `CD3_TILS_COUNTS` trends lower in High `GMSig-4` (median `42.5` vs `60.0`) but is underpowered and non-significant ($\rho = -0.3120$, $p = 0.2994$; Wilcoxon $p = 0.3464$, $n = 13$). Composite `GMSig-4` does **not** separate `ENO1` protein ($\rho = -0.0274$, $p = 0.8901$) or the disjoint 4-protein glycolysis proxy (`LDHA/PKM/GAPDH/PGK1`: $\rho = -0.0049$, $p = 0.9802$), because `GMSig-4` variance in CPTAC TNBC is dominated by the 3-gene immune arm (`SD = 0.4305` vs `0.1849` for `ENO1` arm), even though `ENO1` protein itself tightly tracks the 4-protein glycolysis proxy ($\rho = +0.7614$, $p = 2.54 \times 10^{-6}$).

---

## 2. Data Provenance & Cohort Definition

All four files were downloaded from cBioPortal GitHub LFS (`https://media.githubusercontent.com/media/cBioPortal/datahub/master/public/brca_cptac_2020/...`, total **26,658,918 B = 25.42 MiB**; full manifest in [`cptac_download_manifest.csv`](file:///D:/Bioinformatics/TNBC_Microbiome/02_results/JTM_upgrade/cptac/cptac_download_manifest.csv)):

| File | Bytes | MiB | SHA256 |
|---|---:|---:|---|
| `data_clinical_patient.txt` | 25,254 | 0.0241 | `d9b0dee7ee144105a3fdd37b3573a76558114ec1d8e6fa7262117657517a1275` |
| `data_clinical_sample.txt` | 14,474 | 0.0138 | `f28568f1608c9f63dfeceed0ea301e29ee2dea0efbd88bd826fa2510434896cd` |
| `data_mrna_seq_fpkm.txt` | 17,773,792 | 16.9504 | `00614131229ee33188d7dd20be12d2ac297c052840d93480681852f44e9d0b19` |
| `data_protein_quantification.txt` | 8,845,398 | 8.4356 | `8bef626ce67ff1f70dd229f6a9f88e3bd99284ea38f82ea8eff9348e72bb6a23` |

- **Primary TNBC Cohort (`n = 28`):** `TNBC_UPDATED_CLINICAL_STATUS == "Positive"` (`23` PAM50 Basal, `5` PAM50 Her2; 100% `ERBB2_PROTEOGENOMIC_STATUS == "Negative"`).
- **Strict Clinical IHC Sensitivity Subset (`n = 16`):** `ER_UPDATED_CLINICAL_STATUS == "Negative"`, `PR_CLINICAL_STATUS == "negative"`, and `ERBB2_UPDATED_CLINICAL_STATUS == "Negative"` (excludes 11 cases with `ERBB2_UPDATED_CLINICAL_STATUS = NA` and 1 equivocal case `X05BR045`).
- **RNA Expression Scale Verification:** `data_mrna_seq_fpkm.txt` is **already log2-transformed** (e.g., TNBC `IFNG` spans `-5.246` to `+2.700`, `IDO1` spans `-2.911` to `+7.970`, `ENO1` spans `+7.038` to `+10.261`). Consistent with `frozen_coefficients.json` (`"already transformed gene expression; gene-wise cohort centering and sample SD (n-1), then weighted sum"`), no second `log2(x+1)` transform was applied. Across `n = 28` TNBC, frozen `GMSig4_Score` spans `-0.6270` to `+0.9706` (median `-0.0765`, SD `0.3989`), yielding exact median-split groups (`score > median` = **14 High**, **14 Low**).

---

## 3. Non-Imputed RNA-to-Protein Concordance (TNBC `n = 28` & Strict IHC `n = 16`)

| Gene | MS Protein Detected (TNBC) | TNBC (`n = 28`) Spearman $\rho$ | TNBC Spearman $p$ | TNBC Pearson $r$ ($p$) | Strict IHC (`n = 16`) Pairs | Strict IHC Spearman $\rho$ ($p$) |
|---|---:|---:|---:|---|---:|---|
| **`IFNG`** | **0 / 28 (0.0%)** | — | — | — | 0 | Absent from MS proteome |
| **`CD274`** | **0 / 28 (0.0%)** | — | — | — | 0 | Absent from MS proteome |
| **`ENO1`** | **28 / 28 (100.0%)** | **+0.4625** | **0.0132** | +0.3788 (0.0469) | 16 | +0.2235 (0.4053) |
| **`IDO1`** | **23 / 28 (82.1%)** | **+0.8389** | **$5.70 \times 10^{-7}$** | +0.7970 ($5.32 \times 10^{-6}$) | 14 | +0.8549 ($9.83 \times 10^{-5}$) |
| **`CXCL10`** | **21 / 28 (75.0%)** | **+0.8805** | **$1.40 \times 10^{-7}$** | +0.7910 ($1.97 \times 10^{-5}$) | 12 | +0.8601 ($3.32 \times 10^{-4}$) |
| **`STAT1`** | **28 / 28 (100.0%)** | **+0.8714** | **$1.60 \times 10^{-9}$** | +0.8404 ($2.21 \times 10^{-8}$) | 16 | +0.8882 ($4.49 \times 10^{-6}$) |
| `LDHA` | 28 / 28 (100.0%) | +0.7564 | $3.21 \times 10^{-6}$ | +0.7028 ($3.04 \times 10^{-5}$) | 16 | +0.6912 (0.0030) |
| `PKM` | 28 / 28 (100.0%) | +0.6256 | $3.70 \times 10^{-4}$ | +0.4979 (0.0070) | 16 | +0.5941 (0.0152) |
| `GAPDH` | 28 / 28 (100.0%) | +0.5517 | 0.0023 | +0.5355 (0.0033) | 16 | +0.6206 (0.0103) |
| `PGK1` | 28 / 28 (100.0%) | +0.6606 | $1.30 \times 10^{-4}$ | +0.7187 ($1.65 \times 10^{-5}$) | 16 | +0.6176 (0.0108) |
| `GZMA` | 28 / 28 (100.0%) | +0.7690 | $1.74 \times 10^{-6}$ | +0.7644 ($2.19 \times 10^{-6}$) | 16 | +0.7206 (0.0016) |
| `PRF1` | 22 / 28 (78.6%) | +0.6443 | 0.0012 | +0.7255 ($1.33 \times 10^{-4}$) | 13 | +0.5989 (0.0306) |
| `CD8A` | 27 / 28 (96.4%)* | +0.7108 | $4.72 \times 10^{-5}$ | +0.7696 ($4.32 \times 10^{-6}$) | 15 | +0.7679 ($8.26 \times 10^{-4}$) |

*\*`CD8A` has 27/28 detected in RNA and 27/28 in protein, yielding 26/28 complete RNA–protein pairs.*

- **MNAR Missingness Evidence:** Comparing RNA levels between protein-detected vs protein-NA TNBC cases shows significantly lower mRNA in protein-NA tumors for `IDO1` (`23` vs `5`: median `4.583` vs `0.185`, Mann–Whitney $p = 0.0128$) and `PRF1` (`22` vs `6`: median `1.291` vs `-0.422`, $p = 0.0283$), with a consistent trend for `CXCL10` (`21` vs `7`: median `6.557` vs `5.761`, $p = 0.0623$).

---

## 4. Frozen RNA `GMSig-4` vs Proteomics, ESTIMATE, and CD3 TILs

### 4.1 Primary TNBC Cohort (`n = 28`; 14 High vs 14 Low)

| Domain | Feature | $n$ (High / Low) | Spearman $\rho$ | Spearman $p$ | Median High | Median Low | Wilcoxon $p$ |
|---|---|---:|---:|---:|---:|---:|---:|
| **GMSig-4 Proteins** | `ENO1` protein | 28 (14 / 14) | -0.0274 | 0.8901 | +1.091 | +1.145 | 0.8362 |
| **GMSig-4 Proteins** | **`IDO1` protein** | **23 (11 / 12)** | **-0.7253** | **$9.01 \times 10^{-5}$** | **-1.371** | **+2.351** | **0.00284** |
| **GMSig-4 Proteins** | **`CXCL10` protein** | **21 (11 / 10)** | **-0.7740** | **$3.84 \times 10^{-5}$** | **-0.570** | **+1.938** | **0.00277** |
| **Clinical / TME** | **`ESTIMATE_IMMUNE_SCORE`** | **28 (14 / 14)** | **-0.7997** | **$3.30 \times 10^{-7}$** | **933.9** | **1884.5** | **$4.78 \times 10^{-5}$** |
| **Clinical / TME** | `ESTIMATE_STROMAL_SCORE` | 28 (14 / 14) | +0.0515 | 0.7949 | +80.9 | +112.1 | 1.0000 |
| **Clinical / TME** | `ESTIMATE_TUMORPURITY` | 28 (14 / 14) | +0.5041 | 0.00623 | 0.728 | 0.628 | 0.00824 |
| **Clinical / TME** | `CIBERSORT_ABSOLUTE_SCORE` | 28 (14 / 14) | -0.8183 | $1.04 \times 10^{-7}$ | 0.634 | 1.377 | $1.25 \times 10^{-4}$ |
| **Clinical / TME** | `XCELL_IMMUNE_SCORE` | 28 (14 / 14) | -0.7860 | $7.18 \times 10^{-7}$ | 0.048 | 0.146 | $1.74 \times 10^{-5}$ |
| **Pathology TILs** | **`CD3_TILS_COUNTS`** | **13 (6 / 7)** | **-0.3120** | **0.2994** | **42.5** | **60.0** | **0.3464** |
| **Glycolysis Proteins** | `LDHA` protein | 28 (14 / 14) | +0.0826 | 0.6759 | +0.544 | +0.131 | 0.5351 |
| **Glycolysis Proteins** | `PKM` protein | 28 (14 / 14) | -0.1385 | 0.4822 | -0.089 | +0.290 | 0.4763 |
| **Glycolysis Proteins** | `GAPDH` protein | 28 (14 / 14) | +0.0219 | 0.9119 | +0.774 | +0.431 | 0.8722 |
| **Glycolysis Proteins** | `PGK1` protein | 28 (14 / 14) | -0.0104 | 0.9581 | +0.173 | +0.506 | 0.8005 |
| **Glycolysis Proteins** | `Glycolysis_4Prot_Proxy` (no `ENO1`) | 28 (14 / 14) | -0.0049 | 0.9802 | +0.197 | +0.085 | 0.8362 |
| **Immune Proteins** | **`STAT1` protein** | **28 (14 / 14)** | **-0.8566** | **$6.03 \times 10^{-9}$** | **-1.127** | **+1.220** | **$2.59 \times 10^{-4}$** |
| **Immune Proteins** | **`GZMA` protein** | **28 (14 / 14)** | **-0.7192** | **$1.61 \times 10^{-5}$** | **-0.904** | **+0.441** | **$4.40 \times 10^{-4}$** |
| **Immune Proteins** | **`PRF1` protein** | **22 (10 / 12)** | **-0.6047** | **0.00287** | **-0.858** | **+0.277** | **0.0443** |
| **Immune Proteins** | **`CD8A` protein** | **27 (13 / 14)** | **-0.6490** | **$2.50 \times 10^{-4}$** | **-1.813** | **+0.010** | **0.00611** |
| **Immune Proteins** | `Immune_2Prot_Complete` (`STAT1+GZMA`) | 28 (14 / 14) | -0.8440 | $1.67 \times 10^{-8}$ | -0.652 | +0.766 | $1.81 \times 10^{-4}$ |
| **Immune Proteins** | `Immune_4Prot_Avail` (`STAT1+GZMA+PRF1+CD8A`) | 28 (14 / 14) | -0.8626 | $3.58 \times 10^{-9}$ | -0.579 | +0.500 | $1.50 \times 10^{-4}$ |
| **Immune Proteins** | `Immune_4Prot_Complete` (all 4 detected) | 22 (10 / 12) | -0.8724 | $1.21 \times 10^{-7}$ | -0.511 | +0.538 | 0.00174 |

### 4.2 Strict Clinical IHC TNBC Sensitivity (`n = 16`, Re-Z-Scaled Within `n = 16`)

All key immune proteomic and ESTIMATE associations persist with equal or stronger effect sizes in the strict IHC ER-/PR-/HER2- subset (`n = 16`, 8 High vs 8 Low):
- `IDO1` protein (`n = 14`): $\rho = -0.7451$, $p = 0.00223$; Wilcoxon $p = 0.0127$.
- `CXCL10` protein (`n = 12`): $\rho = -0.9021$, $p = 6.00 \times 10^{-5}$; Wilcoxon $p = 0.00253$.
- `ESTIMATE_IMMUNE_SCORE` (`n = 16`): $\rho = -0.8412$, $p = 4.44 \times 10^{-5}$; Wilcoxon $p = 0.00187$.
- `STAT1` protein (`n = 16`): $\rho = -0.8618$, $p = 1.78 \times 10^{-5}$; Wilcoxon $p = 0.00109$.
- `GZMA` protein (`n = 16`): $\rho = -0.7529$, $p = 7.61 \times 10^{-4}$; Wilcoxon $p = 3.11 \times 10^{-4}$.
- `Immune_4Prot_Avail` (`n = 16`): $\rho = -0.9382$, $p = 7.80 \times 10^{-8}$; Wilcoxon $p = 1.55 \times 10^{-4}$.
- `ENO1` protein (`n = 16`): $\rho = +0.1824$, $p = 0.4991$ (Wilcoxon $p = 0.3823$); `Glycolysis_4Prot_Proxy`: $\rho = +0.2882$, $p = 0.2790$ (Wilcoxon $p = 0.5054$); `CD3_TILS_COUNTS` (`n = 11`): $\rho = -0.2294$, $p = 0.4975$ (Wilcoxon $p = 0.4455$).

---

## 5. Biological & Methodological Interpretation for JTM Manuscript

1. **What CPTAC Validated:** In an independent prospective TNBC cohort (`n = 28`, 0% TCGA overlap), frozen RNA `GMSig-4` exhibits strong, reproducible protein-level concordance with the IFN-$\gamma$/cytotoxic immune axis—not only for its own detected immune genes (`IDO1`, `CXCL10`) but also for completely disjoint cytotoxic/IFN effector proteins (`STAT1`, `GZMA`, `CD8A`, `PRF1`) and bulk immune deconvolution (`ESTIMATE`, `CIBERSORT`, `xCell`), while remaining orthogonal to bulk stromal content (`ESTIMATE_STROMAL_SCORE`, $\rho = 0.0515$, $p = 0.795$).
2. **Why `ENO1` / Glycolysis Proteins Were Null Against Composite `GMSig-4`:** `ENO1` protein strongly tracks core glycolytic enzymes (`LDHA/PKM/GAPDH/PGK1`, $\rho = +0.7614$, $p = 2.54 \times 10^{-6}$) and shows modest positive mRNA–protein coupling ($\rho = +0.4625$, $p = 0.0132$). However, within this small `n = 28` cohort, composite `GMSig-4` variance is overwhelmingly driven by the three co-expressed immune genes (`GMSig4_Immune3_Arm` SD `0.4305`, $\rho = +0.9130$ with `GMSig-4`) relative to `ENO1` (`GMSig4_ENO1_Arm` SD `0.1849`, $\rho = +0.0624$ with `GMSig-4`), and `RNA_ENO1` is weakly positively correlated with `RNA_IFNG` ($\rho = +0.2863$, $p = 0.140$). Thus, composite `GMSig-4` in CPTAC-BRCA should be presented strictly as orthogonal protein-level validation of its immune-effector axis and individual RNA–protein fidelity, **not** as proteomic validation of a 4-protein signature, glycolysis-stratification, or clinical survival.
