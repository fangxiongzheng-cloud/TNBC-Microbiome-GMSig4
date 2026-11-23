# GMSig-4 in triple-negative breast cancer

**A Frozen Four-Gene Tumor Expression Score (GMSig-4) Derived from Microbiome-Linked Host Targets Predicts Survival and Pan-Compartment Immune Depletion in Triple-Negative Breast Cancer**

Fangxiong Zheng, Jinxi Wang, Jiayi Li, Tengfei Ma, Yuchen Zhang, Shuangzhe Li, and Zhongxin Lu

This repository accompanies the manuscript prepared for the *Journal of Translational Medicine*. It contains the analysis and figure scripts for the completed study and the summary tables used to report the results. The workflow integrates gut microbiome profiles, germline genetic evidence, bulk tumor expression and clinical outcomes, and single-cell and structural analyses.

The four-gene Ridge Cox score reported in the study is

```text
GMSig-4 = -0.2518*Z_IFNG - 0.1509*Z_IDO1 - 0.0719*Z_CXCL10 + 0.1849*Z_ENO1
```

Each `Z` is gene expression standardized within the cohort being scored. The displayed coefficients are rounded to four decimals. The analysis script fits the model in TCGA-TNBC using the original five-fold cross-validation and seed of 2026, then reuses the fitted coefficients at full precision for validation. The development model used 29 events for four predictors, yielding 7.25 events per variable (EPV).

## JTM release v1.1.0

The [v1.1.0 release](https://github.com/fangxiongzheng-cloud/TNBC-Microbiome-GMSig4/releases/tag/v1.1.0) archives the scripts, frozen model and processed results used in the JTM revision.

| Analysis | Where to look |
| --- | --- |
| Frozen model | [Full-precision coefficients](02_results/JTM_upgrade/model_freeze/frozen_coefficients.json), [scoring function](00_scripts/26a_gmsig4_score.R) and [SHA-256 manifest](02_results/JTM_upgrade/model_freeze/sha256_manifest.csv). The manifest records the original analysis workspace. Its 16 model-freeze files and three `26a` scripts can be verified byte-for-byte in this release. Its other entries describe seven historical scripts that were adapted for the public repository and seven model or score files kept in the analysis workspace. |
| External validation | `26b` handles eligibility and harmonization, `26c` evaluates the frozen score, and `26i` plots the results. See [eligibility](02_results/JTM_upgrade/cohort_eligibility) and [validation tables](02_results/JTM_upgrade/independent_validation) for METABRIC, GSE21653 and GSE25066. |
| Spatial, proteomic and docking analyses | `26d`–`26e` audit structures and available data, `26g` and `26h` analyze GSE210616 and CPTAC, and `26j` and `26k` plot the results. See the [spatial tables](02_results/JTM_upgrade/spatial_gse210616), [CPTAC results](02_results/JTM_upgrade/cptac) and [docking audit](02_results/JTM_upgrade/docking). |
| Reconciliation and method audit | `26f` reconciles the historical tables. The [reconciled tables](02_results/JTM_upgrade/reconciled_tables), [method audit](02_results/JTM_upgrade/method_audit) and [PRJNA383849 subject-level recomputation summary](02_results/PRJNA383849_subject_collapsed_recomputation.csv) record the checks behind the revision. |

For GCST010098, GCST010100 and GCST90454344, the eQTLGen tables pool the primary inverse-variance weighted and Wald ratio tests in one Benjamini–Hochberg adjustment per outcome. The `fdr` column contains the revised values and `fdr_by_method_old` preserves the previous values. Each `MR_eQTLGen_*_results_unifiedBH.csv` is byte-identical to the corresponding `MR_eQTLGen_*_results.csv` used by the analysis scripts. The Figure 2 script applies the same per-outcome adjustment when combining endpoints.

The [Figure 2 PDF](03_figures/Figure2.pdf) and [PNG](03_figures/Figure2.png) are the submission exports. The [plotting script](00_scripts/22_fig2_mr_coloc_mediation.R) reproduces the revised panel D colocalization posteriors and the panel E GTEx FDR annotation (0.121).

The archived `14c`, `14e`, `14h` and several `26` scripts retain paths from the original Windows/WSL analysis environment. Review the path settings before running them elsewhere. Use a working copy when regenerating outputs to preserve the released model and its hashes.

## Local setup

1. Start from the repository root. The portable scripts use `PROJECT_ROOT` and default to the current directory. The archived JTM scripts have the path settings described above. To run from another working directory, set `PROJECT_ROOT` in the shell to the path of your analysis checkout.
2. Use **R 4.4.3** with **Bioconductor 3.20** and **Python 3.12**. Keep `Rscript` and `python` available in the active environment. Use Linux or WSL for docking if your platform lacks compatible Vina wheels. [Bioconductor releases](https://bioconductor.org/about/release-announcements/)
3. Generate your own token in the OpenGWAS profile. Copy `.Renviron.example` to `.Renviron` in the analysis checkout and replace the JWT placeholder locally. Use `PROJECT_ROOT=.` when launching from the repository root and the checkout path when launching from elsewhere. Start a fresh R process after editing. Keep the token and source inputs local. [OpenGWAS authentication guide](https://mrcieu.github.io/ieugwasr/articles/guide.html#authentication)
4. Stage the public inputs under `01_data/` as described below. The released `02_results/` contains summary tables and selected derived score tables from public cohorts. Run in a separate checkout to preserve that snapshot because analyses write new tables there. Figures are generated under `03_figures/`.

### Dependencies

| Runtime | Core packages |
| --- | --- |
| R data and graphics | `data.table`, `jsonlite`, `ggplot2`, `ggrepel`, `patchwork`, `ragg`, `systemfonts`, `scales` |
| R microbiome and genetics | `vegan`, `MicrobiomeStat`, `TwoSampleMR`, `ieugwasr`, `coloc`, and `MRPRESSO` for the MR-PRESSO sensitivity step |
| R survival and classification | `survival`, `glmnet`, `randomForestSRC`, `rms`, `pROC`, `timeROC` |
| R expression and pathways | `GEOquery`, `limma`, `GSVA`, `fgsea`, `clusterProfiler`, `AnnotationDbi`, `hgu133plus2.db`, `org.Hs.eg.db`, `msigdbr` |
| Python data and plots | `numpy`, `pandas`, `scipy`, `scikit-learn`, `requests`, `matplotlib`, `Pillow`, `pyliftover` |
| Python single-cell and docking | `scanpy`, `anndata`, `h5py`, `rdkit`, `meeko`, `gemmi`, and `vina==1.2.7` |

Run `Rscript 00_scripts/00_install_r_packages.R` with the intended R active. `00_run_install_wsl.sh` is a Bash wrapper for the same installer. Install `MRPRESSO` from its upstream source to reproduce that sensitivity step. [MR-PRESSO](https://github.com/rondolab/MR-PRESSO)

```sh
python -m pip install numpy pandas scipy scikit-learn requests matplotlib Pillow pyliftover scanpy anndata h5py rdkit meeko gemmi
python -m pip install vina==1.2.7
```

The curl-based downloader uses `curl` from `PATH`. For docking, expose the Meeko `mk_prepare_receptor.py` command on `PATH` or set `MEEKO_RECEPTOR` to the executable path. The Arial figure layouts can use system fonts or four local Arial faces supplied via `ARIAL_FONT_DIR`. [AutoDock Vina documentation](https://autodock-vina.readthedocs.io/en/latest/docking_requirements.html)

## Execution order

Run Module 1, the first two rows of Module 3, Module 2, the remaining rows of Module 3, and Module 4. In each row run scripts left to right, using `python 00_scripts/<script>.py` or `Rscript 00_scripts/<script>.R`.

### 1. HGMT and GMrepo gut microbiome analysis

| Order | Scripts | Output or prerequisite |
| --- | --- | --- |
| 1 | `01_download_gutmgene_and_gmrepo_meta.py` → `02_explore_gutmgene_summary.py` | Human gutMGene tables and candidate target genes. Stage the GMrepo breast-neoplasm project response as described below. |
| 2 | `03_fetch_gmrepo_abundances.py` → `04_gmrepo_diff_analysis.R` → `05_map_gmrepo_to_gutmgene.py` | QC-passed run profiles, differential abundance and diversity, and GMrepo-supported gutMGene targets. |
| 3 | `15_hgmt_analysis.py` → `16_hgmt_kegg_pathway_enrichment.py` → `15b_hgmt_all8_cohorts.py` | Case-control profiles and KO matrix, KEGG pathway matrix, then the eight-project analysis and the local Figure 1B input. |
| 4 | `22_fig1_microbiome_hgmt_v2.R` | Figure 1 uses locally regenerated PRJNA658160 diversity and abundance data. `22_fig1_microbiome_hgmt.R` draws the summary-table version. |

### 2. Subtype-stratified bidirectional MR, Steiger tests, colocalization, and mediation

| Order | Scripts | Output or prerequisite |
| --- | --- | --- |
| 1 | `06b_fetch_mibiogen_tophits.py` → `14_mibiogen_bcac_mr.R` → `14b_mibiogen_bcac_sensitivity.R` | MiBioGen instruments, forward BCAC MR, heterogeneity, pleiotropy, and sensitivity. |
| 2 | `17_eQTLGen_OpenGWAS_SMR_coloc.R` | The expression evidence matrix from Module 3 supports instrument selection, eQTLGen MR, and colocalization. |
| 3 | `06i_download_tbi_indexes.py` → `06j_liftover_ivs_grch38.py` → `06k_tbi_range_extract.py` | Tabix indexes, GRCh37-to-GRCh38 mapping, and variant extraction for three GWAS Catalog studies. `06j_map_ivs_grch38.py` is an API-based mapping utility. |
| 4 | `14c_gcst010100_tnbc_mr.R` → `14e_gcst90454344_puretnbc_mr.R` → `14h_gcst010098_overallbc_mr.R` → `14f_cross_subtype_mr_master.R` | TNBC, pure TNBC, overall breast cancer, and cross-subtype comparisons. |
| 5 | `06l_get_gene_coords.py` → `06m_extract_gtex_breast_and_tnbc_regions.py` → `14i_gtex_breast_mr_and_tnbc_coloc.R` → `14g_mediation_mr_taxa_to_eqtlgen.R` | GTEx breast and TNBC regional extraction, tissue-specific MR and colocalization, and two-step mediation. If GTEx extraction reports missing records, run `06n_retry_missing_gtex_breast.py` before `14i`. |
| 6 | `25a_lookup_reverse_mr_ids.R` → `25d1_fetch_bc_ivs.R` → `25d2_extract_tnbc_reverse_ivs.py` → `25d3_reverse_mr_steiger_figS1.R` → `22_fig2_mr_coloc_mediation.R` | Reverse MR, forward Steiger directionality tests, Supplementary Figure S1, and the MR and colocalization figure. |

### 3. GMSig-4 survival, immunotherapy response, and decision curves

| Order | Scripts | Output or prerequisite |
| --- | --- | --- |
| 1 | `06_download_core_datasets.py` → `06d_range_downloader.py` → `06f_download_gse96058_scrna.py` | TCGA/Xena and GEO inputs. `06c_download_transcriptomics.py` supports curl-based retries for applicable files. |
| 2 | `07_tcga_tnbc_gutmgene_screen.R` → `09_ispy2_pembrolizumab_screen.R` → `11_validate_gse58812_gse103091_gse173839.R` | Multi-cohort screening and evidence matrix used by Module 2. |
| 3 | `18_signature_candidate_audit.py` → `19a_prepare_cohort_matrices.R` → `19a2_prepare_gse96058.R` | Integrated candidate evidence and locally cached cohort matrices. |
| 4 | `19b_ml_compact_signature.R` → `19c_screen_and_benchmark_all5.R` → `19d_benchmark_all5_models.R` → `19e_compare_signature_candidates.R` | Compact signatures and survival and immunotherapy benchmarks across five cohorts. |
| 5 | `19f_signature_clinical_nomogram.R` → `19g_tme_immune_ici_analysis.R` → `19h_calibration_and_nomogram_points.R` → `25b_clinical_dca.R` | GMSig-4 fitting and survival validation, pCR and tumor microenvironment analyses, calibration, and decision curve analysis (DCA). This row regenerates model objects and patient scores locally. |
| 6 | `22_fig3_signature_prognosis.R` → `22_fig4_tme_gsea.R` → `24_preflight_checks.R` | Prognosis, calibration, and pathway figures, followed by the subgroup scoring needed for Figure 5. |

### 4. Single-cell communication and molecular docking

| Order | Scripts | Output or prerequisite |
| --- | --- | --- |
| 1 | `20a_scrna_tnbc_atlas.py` → `20b_scrna_tnbc_summary_and_tme.py` → `25c_scrna_cell_comm.py` | The GSE176078 TNBC subset of 42,512 cells from 10 tumors, cell-type and pseudobulk summaries, and permutation-based ligand–receptor communication. Extract the atlas first. |
| 2 | `22_fig5_ici_and_scrna.py` | Figure 5 uses the ICI results from `19g` and `24`, the UMAP from `20a`, and the cell-type summaries from `20b`. |
| 3 | `21a_fetch_structures_and_ligands.py` → `21b_autodock_vina_pipeline.py` → `22_fig6_docking_and_mechanism.py` | Structure and ligand manifests, AutoDock Vina 1.2.7 affinities and poses, contact residues, and the data-driven Figure 6. |
| 4 | `22_fig6_docking_and_mechanism_v2.py` (optional) | The illustrated layout requires local pocket renders named `<gene>_<metabolite>.png` in `01_data/docking/figure_assets/pockets3d/` and `panel6d_background.jpg` in `figure_assets/`. |

The `00_scripts/audit_release.py` checker describes the two-directory v1.0.0 archive. Its allowlist predates the JTM model files, figures and nested result directories. After recomputation, use `23_audit_all_results.py` to inspect the generated analysis outputs.

## Public data accessions

| Source | Accessions or release | Use |
| --- | --- | --- |
| HGMT, 16S | [PRJNA345373](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA345373), [PRJNA383849](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA383849), [PRJNA566060](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA566060), [PRJNA658160](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA658160), [PRJNA751260](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA751260) | Taxonomic profiles and phenotype strata. [HGMT paper and data availability](https://pmc.ncbi.nlm.nih.gov/articles/PMC12642048/) |
| HGMT, shotgun WGS | [PRJNA453965](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA453965), [PRJNA494824](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA494824), [PRJNA718520](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA718520) | Taxonomic and KEGG orthology profiles. [HGMT portal](https://mai.fudan.edu.cn/hgmt) |
| GMrepo | [GMrepo v2](https://pmc.ncbi.nlm.nih.gov/articles/PMC8728112/), projects PRJNA383849, PRJNA658160, and [PRJNA872152](https://www.ncbi.nlm.nih.gov/bioproject/PRJNA872152) | Breast Neoplasms (D001943) and Health (D006262) profiles. [GMrepo portal](https://gmrepo.humangut.info/) |
| gutMGene | [gutMGene v2.0](https://bio-computing.hrbmu.edu.cn/gutMGene/) | Human microbe–gene, microbe–metabolite, and metabolite–gene associations. |
| MiBioGen | [MiBioGen top hits](https://mibiogen.gcc.rug.nl/) | Microbial trait GWAS instruments. |
| OpenGWAS / BCAC | [ieu-a-1126](https://opengwas.io/datasets/ieu-a-1126), [ieu-a-1128](https://opengwas.io/datasets/ieu-a-1128), [ieu-a-1135](https://opengwas.io/datasets/ieu-a-1135) | Overall breast cancer, pooled ER-negative breast cancer, and OncoArray ER-negative breast cancer, respectively. |
| GWAS Catalog | [GCST010098](https://www.ebi.ac.uk/gwas/studies/GCST010098), [GCST010100](https://www.ebi.ac.uk/gwas/studies/GCST010100), [GCST90454344](https://www.ebi.ac.uk/gwas/studies/GCST90454344) | Overall breast cancer, TNBC, and the pure-TNBC endpoint used in the analysis, respectively. |
| eQTLGen | [eQTLGen cis-eQTL](https://www.eqtlgen.org/cis-eqtls.html) | Blood cis-eQTL instruments and OpenGWAS `eqtl-a-ENSG…` gene datasets. |
| GTEx | [GTEx v8](https://gtexportal.org/home/downloads/adult-gtex/bulk_tissue_expression), Breast Mammary Tissue | Tissue-specific cis-eQTLs via indexed eQTL Catalogue files. [eQTL Catalogue downloads](https://www.ebi.ac.uk/eqtl/Data_access/) |
| TCGA | [TCGA-BRCA](https://portal.gdc.cancer.gov/projects/TCGA-BRCA) | Development and internal survival analysis using UCSC Xena HiSeqV2, phenotype, and survival files. |
| GEO survival | [GSE58812](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE58812), [GSE96058](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE96058) | External survival validation, including SCAN-B. |
| GEO immunotherapy | [GSE194040](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE194040), [GSE173839](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE173839) | I-SPY2 pembrolizumab and durvalumab/olaparib pCR analyses. |
| GEO single-cell | [GSE176078](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE176078) | Wu et al. breast cancer atlas. Use the processed GEO expression matrix and metadata. |
| GEO exploratory screen | [GSE103091](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE103091) | Expression cohort inspection. |

### Local input layout

- **HGMT.** Obtain the project details, run annotations, bacterial profiles, and WGS KO profiles through the HGMT data downloads. Under `01_data/HGMT/`, save `<PRJNA>_details.json` as an object, `<PRJNA>_runs.json` as an array of run annotations, and `<PRJNA>_feature_bacteria.json` as an array with `run_id`, `taxa`, and `abundance`. For the three WGS studies, also supply `<PRJNA>_feature_ko.json` with `run_id`, `ko`, and `abundance`. Preserve phenotype and menopause annotations. Save KEGG pathway links, pathway names, and KO names as `kegg_ko_pathway_link.tsv`, `kegg_pathway_names.tsv`, and `kegg_ko_names.tsv`. [KEGG REST API](https://www.kegg.jp/kegg/rest/keggapi.html)
- **GMrepo.** Save the raw `getAssociatedProjectsByMeshID` response for `D001943` as `01_data/GMrepo/getAssociatedProjectsByMeshID_D001943.json`. Save each of the three GMrepo projects' run exports as `all_runs_in_project_<PRJNA>.tsv` in that directory. Retain `QCstatus`, `run_id`, `project_id`, `sample_id`, `disease`, and `term`.
- **GWAS and eQTL.** Put the UCSC `hg19ToHg38.over.chain.gz` file and GWAS Catalog indexes in `01_data/GWAS_MR/`. The GTEx v8 breast index and extracted regions go under `01_data/eQTL/`. Scripts `06i`, `06k`, and `06m` retrieve indexes or regions from EBI. If an upstream filename has changed, use the matching harmonized file and index for the same accession and genome build. [UCSC liftOver chains](https://hgdownload.soe.ucsc.edu/goldenPath/hg19/liftOver/)
- **TCGA and GEO.** Keep bulk expression, series matrices, and phenotype or biomarker files under `01_data/TCGA_GEO/`, using the names in the download scripts. Supply the GSE194040 probe annotation as `GSE194040_ProbeAnnotation.txt.gz` for probe-level mapping. For pathway analyses, save the human-symbol MSigDB Hallmark and KEGG gene sets as `01_data/TCGA_GEO/msigdb/hallmark.gmt` and `01_data/TCGA_GEO/msigdb/kegg.gmt`. [MSigDB collections](https://www.gsea-msigdb.org/gsea/msigdb/human/collections.jsp)
- **Single-cell and docking.** Extract the GSE176078 processed GEO supplement into `01_data/scRNA/Wu_etal_2021_BRCA_scRNASeq/` with `metadata.csv`, `count_matrix_genes.tsv`, `count_matrix_barcodes.tsv`, and `count_matrix_sparse.mtx` at that level. Script `21a` stages structures and ligands under `01_data/docking/`. The released manifests record the structures and compounds used in the study.

## Citation and contact

Use the repository citation for the full title and author order. [CITATION.cff](CITATION.cff)

The corresponding author is Zhongxin Lu, Department of Medical Laboratory, The Central Hospital of Wuhan, Tongji Medical College, Huazhong University of Science and Technology, Wuhan 430014, Hubei, China. [luzhongxin@zxhospital.com](mailto:luzhongxin@zxhospital.com)

The code is available under the MIT License. Third-party inputs are governed by their source terms. [LICENSE](LICENSE)
