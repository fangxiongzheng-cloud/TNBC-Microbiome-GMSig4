#!/usr/bin/env Rscript
# Reconcile historical tables per Stage 0/A method_audit without mutating original provenance files.
# Writes clean audited tables to 02_results/JTM_upgrade/reconciled_tables/.

.libPaths(c("/root/miniconda3/envs/bioinfo-r-core/lib/R/library", "/root/R/tnbc_microbiome_lib", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(pROC)
  library(jsonlite)
})

root    <- "/mnt/d/Bioinformatics/TNBC_Microbiome"
res_dir <- file.path(root, "02_results")
cache   <- file.path(root, "01_data/TCGA_GEO/clean_cache")
out_dir <- file.path(res_dir, "JTM_upgrade/reconciled_tables")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 1. Reconcile 19f Prognostic Summary: true per-SD HR vs raw-unit HR
cat("=== 1. Reconciling 19f survival per-SD vs raw-unit HRs ===\n")
source(file.path(root, "00_scripts/26a_gmsig4_score.R"))
m <- readRDS(file.path(res_dir, "GMSig4_Model_Objects.rds"))
w_ridge <- m$weights

tcga <- readRDS(file.path(cache, "TCGA_BRCA_clean.rds"))
g58  <- readRDS(file.path(cache, "GSE58812_clean.rds"))
g960 <- readRDS(file.path(cache, "GSE96058_TNBC_clean.rds"))
g173 <- readRDS(file.path(cache, "GSE173839_clean.rds"))
g194 <- readRDS(file.path(cache, "GSE194040_clean.rds"))

eval_reconciled <- function(cname, time_yr, event, raw_score) {
  y <- Surv(time_yr, event)
  sd_val <- sd(raw_score)
  z_score <- as.numeric(scale(raw_score))
  cx_raw <- summary(coxph(y ~ raw_score))
  cx_sd  <- summary(coxph(y ~ z_score))
  grp <- factor(ifelse(raw_score > median(raw_score), "High", "Low"), levels = c("Low", "High"))
  cx_grp <- summary(coxph(y ~ grp))
  km_p <- 1 - pchisq(survdiff(y ~ grp)$chisq, 1)
  # Reverse KM median follow-up
  rkm <- survfit(Surv(time_yr, 1 - event) ~ 1)
  rev_km_med <- unname(summary(rkm)$table["median"])
  data.table(
    cohort = cname, n = length(time_yr), events = sum(event),
    obs_median_time_yr = median(time_yr),
    reverse_KM_followup_yr = rev_km_med,
    score_sd = sd_val,
    C_index = concordance(y ~ raw_score, reverse = TRUE)$concordance,
    KM_logrank_p = km_p,
    HighLow_HR = cx_grp$conf.int[1, "exp(coef)"],
    HighLow_LCI = cx_grp$conf.int[1, "lower .95"],
    HighLow_UCI = cx_grp$conf.int[1, "upper .95"],
    HR_per_raw_unit = cx_raw$conf.int[1, "exp(coef)"],
    HR_per_raw_LCI  = cx_raw$conf.int[1, "lower .95"],
    HR_per_raw_UCI  = cx_raw$conf.int[1, "upper .95"],
    HR_per_cohort_SD = cx_sd$conf.int[1, "exp(coef)"],
    HR_per_SD_LCI    = cx_sd$conf.int[1, "lower .95"],
    HR_per_SD_UCI    = cx_sd$conf.int[1, "upper .95"],
    Continuous_Cox_p = cx_sd$coefficients[1, "Pr(>|z|)"]
  )
}

idx_t  <- which(tcga$clin$is_TNBC & !is.na(tcga$clin$OS) & !is.na(tcga$clin$OS.time) & tcga$clin$OS.time > 0)
g58_i  <- which(!is.na(g58$clin$os_event) & !is.na(g58$clin$os_days) & g58$clin$os_days > 0)
g960_s <- which(g960$clin$is_TNBC_strict & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_e <- which(g960$clin$is_TNBC_ext    & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_b <- which(g960$clin$pam50_subtype == "Basal" & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))

s_t    <- score_gmsig4(tcga$expr, w_ridge, colnames(tcga$expr)[idx_t])$GMSig4
s_58   <- score_gmsig4(g58$expr,  w_ridge, colnames(g58$expr)[g58_i])$GMSig4
s_960s <- score_gmsig4(g960$expr, w_ridge, colnames(g960$expr)[g960_s])$GMSig4
s_960e <- score_gmsig4(g960$expr, w_ridge, colnames(g960$expr)[g960_e])$GMSig4
s_960b <- score_gmsig4(g960$expr, w_ridge, colnames(g960$expr)[g960_b])$GMSig4

surv_rec <- rbindlist(list(
  eval_reconciled("TCGA_TNBC_OS",            tcga$clin$OS.time[idx_t]/365.25,  tcga$clin$OS[idx_t],          s_t),
  eval_reconciled("GSE58812_TNBC_OS",        g58$clin$os_days[g58_i]/365.25,   g58$clin$os_event[g58_i],     s_58),
  eval_reconciled("GSE58812_TNBC_MFS",       g58$clin$mfs_days[g58_i]/365.25,  g58$clin$meta_event[g58_i],   s_58),
  eval_reconciled("GSE96058_TNBC_Strict_OS", g960$clin$os_days[g960_s]/365.25, g960$clin$os_event[g960_s],   s_960s),
  eval_reconciled("GSE96058_TNBC_Ext_OS",    g960$clin$os_days[g960_e]/365.25, g960$clin$os_event[g960_e],   s_960e),
  eval_reconciled("GSE96058_Basal_OS",       g960$clin$os_days[g960_b]/365.25, g960$clin$os_event[g960_b],   s_960b)
))
fwrite(surv_rec, file.path(out_dir, "Reconciled_Prior_Survival_Summary.csv"))
print(surv_rec[, .(cohort, n, events, HR_per_raw_unit = round(HR_per_raw_unit, 3), HR_per_cohort_SD = round(HR_per_cohort_SD, 3), reverse_KM_followup_yr = round(reverse_KM_followup_yr, 2))])

# Multivariable per-SD vs raw-unit check on TCGA & GSE96058_Ext
tcga_df <- fread(file.path(res_dir, "GMSig4_PatientScores_TCGA_TNBC.csv"))
tcga_cc <- tcga_df[complete.cases(tcga_df[, .(OS_yr, OS, GMSig4, Age, Stage_III_IV, Node_Pos)])]
tcga_cc[, GMSig4_z := as.numeric(scale(GMSig4))]
g960_df <- fread(file.path(res_dir, "GMSig4_PatientScores_GSE96058_TNBC_Ext.csv"))
g960_cc <- g960_df[complete.cases(g960_df[, .(OS_yr, OS, GMSig4, Age, Tumor_mm, Node_Pos)])]
g960_cc[, GMSig4_z := as.numeric(scale(GMSig4))]

mv_extract <- function(fit, coh, unit_lbl) {
  s <- summary(fit)
  data.table(cohort = coh, gmsig4_unit = unit_lbl, variable = rownames(s$coefficients),
             n = fit$n, events = fit$nevent,
             HR = s$conf.int[, "exp(coef)"], LCI = s$conf.int[, "lower .95"],
             UCI = s$conf.int[, "upper .95"], pval = s$coefficients[, "Pr(>|z|)"])
}
mv_rec <- rbindlist(list(
  mv_extract(coxph(Surv(OS_yr, OS) ~ GMSig4   + Age + Stage_III_IV + Node_Pos, data = tcga_cc), "TCGA_TNBC", "per_raw_unit"),
  mv_extract(coxph(Surv(OS_yr, OS) ~ GMSig4_z + Age + Stage_III_IV + Node_Pos, data = tcga_cc), "TCGA_TNBC", "per_cohort_SD"),
  mv_extract(coxph(Surv(OS_yr, OS) ~ GMSig4   + Age + Tumor_mm     + Node_Pos, data = g960_cc), "GSE96058_TNBC_Ext", "per_raw_unit"),
  mv_extract(coxph(Surv(OS_yr, OS) ~ GMSig4_z + Age + Tumor_mm     + Node_Pos, data = g960_cc), "GSE96058_TNBC_Ext", "per_cohort_SD")
))
fwrite(mv_rec, file.path(out_dir, "Reconciled_Prior_Multivariable_Cox.csv"))

# 2. Reconcile GSE173839 single-gene STAT1/CD274/PDCD1 vs trial STAT1_sig/PDL1/PD1
cat("\n=== 2. Reconciling GSE173839 single-gene vs trial biomarker AUCs ===\n")
auc_safe <- function(y, x) as.numeric(auc(roc(y, x, direction = "<", quiet = TRUE)))
hr_col <- grep("^HR\\.status", colnames(g173$clin), value = TRUE)[1]
arms173 <- list(
  GSE173839_DurvaOlap_All  = which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1)),
  GSE173839_DurvaOlap_TNBC = which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0),
  GSE173839_Ctrl_All       = which(g173$clin$Arm == "control"             & g173$clin$pCR.status %in% c(0, 1)),
  GSE173839_Ctrl_TNBC      = which(g173$clin$Arm == "control"             & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0)
)
ici173_rec <- rbindlist(lapply(names(arms173), function(nm) {
  ii <- arms173[[nm]]
  y  <- g173$clin$pCR.status[ii]
  sc <- score_gmsig4(g173$expr, w_ridge, colnames(g173$expr)[ii])$GMSig4
  data.table(
    cohort = nm, n = length(ii), pCR_n = sum(y),
    GMSig4_AUC = auc_safe(y, -sc),
    SingleGene_CD274_AUC = auc_safe(y, g173$expr["CD274", ii]),
    SingleGene_PDCD1_AUC = auc_safe(y, g173$expr["PDCD1", ii]),
    SingleGene_STAT1_AUC = auc_safe(y, g173$expr["STAT1", ii]),
    TrialBiomarker_PDL1_AUC     = auc_safe(y, as.numeric(g173$clin$PDL1[ii])),
    TrialBiomarker_PD1_AUC      = auc_safe(y, as.numeric(g173$clin$PD1[ii])),
    TrialBiomarker_STAT1sig_AUC = auc_safe(y, as.numeric(g173$clin$STAT1_sig[ii]))
  )
}))
fwrite(ici173_rec, file.path(out_dir, "Reconciled_GSE173839_SingleGene_vs_TrialBiomarkers.csv"))
print(ici173_rec)

# 3. Reconcile 14i GTEx Breast MR FDR across all 19 primary tests + 17 eQTLGen coloc summary
cat("\n=== 3. Reconciling GTEx MR FDR (all 19 primary genes) & eQTLGen coloc ===\n")
gt_mr <- fread(file.path(res_dir, "MR_GTExBreast_GCST010100_TNBC_results.csv"))
gt_prim <- gt_mr[method %in% c("Inverse variance weighted", "Wald ratio")]
gt_prim[, fdr_by_method_old := p.adjust(pval, method = "BH"), by = method]
gt_prim[, fdr_all19_genes   := p.adjust(pval, method = "BH")]
fwrite(gt_prim, file.path(out_dir, "Reconciled_GTExBreast_MR_All19_FDR.csv"))
print(gt_prim[order(pval)][1:5, .(exposure, method, nsnp, OR = round(OR, 3), pval = signif(pval, 4), fdr_by_method_old = round(fdr_by_method_old, 4), fdr_all19_genes = round(fdr_all19_genes, 4))])

eq_coloc <- fread(file.path(res_dir, "Coloc_eQTLGen_BCAC_results.csv"))
eq_coloc[, PP.H3_plus_H4 := PP.H3 + PP.H4]
eq_coloc[, cond_H4_heuristic_ge050 := ifelse(PP.H3_plus_H4 >= 0.50, PP.H4 / PP.H3_plus_H4, NA_real_)]
fwrite(eq_coloc, file.path(out_dir, "Reconciled_eQTLGen_Coloc_Summary.csv"))
cat("eQTLGen coloc rows with PP.H3+PP.H4 >= 0.50:", sum(eq_coloc$PP.H3_plus_H4 >= 0.50), "/", nrow(eq_coloc), "\n")
cat("DONE.\n")
