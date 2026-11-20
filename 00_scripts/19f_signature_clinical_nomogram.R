project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(timeROC)
  library(pROC)
})

cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
res_dir   <- file.path(project_root, "02_results")

tcga <- readRDS(file.path(cache_dir, "TCGA_BRCA_clean.rds"))
g58  <- readRDS(file.path(cache_dir, "GSE58812_clean.rds"))
g960 <- readRDS(file.path(cache_dir, "GSE96058_TNBC_clean.rds"))
g194 <- readRDS(file.path(cache_dir, "GSE194040_clean.rds"))
g173 <- readRDS(file.path(cache_dir, "GSE173839_clean.rds"))

# Locked primary 4-gene signature (GMSig-4): IFNG + IDO1 + CXCL10 + ENO1
sig_genes <- c("IFNG", "IDO1", "CXCL10", "ENO1")

idx_t <- which(tcga$clin$is_TNBC & !is.na(tcga$clin$OS) & !is.na(tcga$clin$OS.time) & tcga$clin$OS.time > 0)
y_tcga <- Surv(tcga$clin$OS.time[idx_t] / 365.25, tcga$clin$OS[idx_t])
Z_tcga <- scale(t(tcga$expr[sig_genes, idx_t, drop = FALSE]))

set.seed(2026)
cv_fit <- cv.glmnet(Z_tcga, y_tcga, family = "cox", alpha = 0, nfolds = 5)
w_ridge <- as.numeric(coef(cv_fit, s = "lambda.min"))
names(w_ridge) <- sig_genes
cat("Locked GMSig-4 Ridge Cox coefficients (trained on TCGA-TNBC n=", length(idx_t), "):\n", sep = "")
print(round(w_ridge, 6))

# Also extract univariate weights for sensitivity checks
w_uni <- sapply(sig_genes, function(g) coef(coxph(y_tcga ~ Z_tcga[, g]))[1])
cat("Univariate Cox betas:\n"); print(round(w_uni, 6))

score_cohort <- function(expr, idx, weights = w_ridge) {
  Z <- scale(t(expr[names(weights), idx, drop = FALSE]))
  as.numeric(Z %*% weights)
}

# 1. Evaluate all survival metrics + 1/3/5-yr timeROC across TCGA, GSE58812 OS/MFS, GSE96058 Strict/Ext/Basal
eval_surv_cohort <- function(cname, time_yr, event, risk_score) {
  y <- Surv(time_yr, event)
  grp <- factor(ifelse(risk_score > median(risk_score), "High", "Low"), levels = c("Low", "High"))
  sd_fit <- survdiff(y ~ grp)
  km_p <- 1 - pchisq(sd_fit$chisq, 1)
  cx <- coxph(y ~ risk_score)
  cx_sum <- summary(cx)
  cx_grp <- summary(coxph(y ~ grp))
  cind <- concordance(y ~ risk_score, reverse = TRUE)$concordance
  troc <- tryCatch(timeROC(T = time_yr, delta = event, marker = risk_score, cause = 1,
                           weighting = "marginal", times = c(1, 3, 5), iid = FALSE), error = function(e) NULL)
  auc1 <- if (!is.null(troc)) troc$AUC[which.min(abs(troc$times - 1))] else NA_real_
  auc3 <- if (!is.null(troc)) troc$AUC[which.min(abs(troc$times - 3))] else NA_real_
  auc5 <- if (!is.null(troc)) troc$AUC[which.min(abs(troc$times - 5))] else NA_real_
  data.table(
    cohort = cname, n = length(time_yr), events = sum(event),
    median_followup_yr = median(time_yr),
    C_index = cind, KM_logrank_p = km_p,
    HighLow_HR = cx_grp$conf.int[1, "exp(coef)"],
    HighLow_LCI = cx_grp$conf.int[1, "lower .95"],
    HighLow_UCI = cx_grp$conf.int[1, "upper .95"],
    HighLow_Cox_p = cx_grp$coefficients[1, "Pr(>|z|)"],
    Continuous_HR_perSD = cx_sum$conf.int[1, "exp(coef)"],
    Continuous_Cox_p = cx_sum$coefficients[1, "Pr(>|z|)"],
    AUC_1yr = auc1, AUC_3yr = auc3, AUC_5yr = auc5
  )
}

s_tcga <- score_cohort(tcga$expr, idx_t)

g58_idx <- which(!is.na(g58$clin$os_event) & !is.na(g58$clin$os_days) & g58$clin$os_days > 0)
s_g58   <- score_cohort(g58$expr, g58_idx)

g960_s <- which(g960$clin$is_TNBC_strict & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_e <- which(g960$clin$is_TNBC_ext    & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_b <- which(g960$clin$pam50_subtype == "Basal" & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
s_960s <- score_cohort(g960$expr, g960_s)
s_960e <- score_cohort(g960$expr, g960_e)
s_960b <- score_cohort(g960$expr, g960_b)

surv_summary <- rbindlist(list(
  eval_surv_cohort("TCGA_TNBC_OS", tcga$clin$OS.time[idx_t] / 365.25, tcga$clin$OS[idx_t], s_tcga),
  eval_surv_cohort("GSE58812_TNBC_OS", g58$clin$os_days[g58_idx] / 365.25, g58$clin$os_event[g58_idx], s_g58),
  eval_surv_cohort("GSE58812_TNBC_MFS", g58$clin$mfs_days[g58_idx] / 365.25, g58$clin$meta_event[g58_idx], s_g58),
  eval_surv_cohort("GSE96058_TNBC_Strict_OS", g960$clin$os_days[g960_s] / 365.25, g960$clin$os_event[g960_s], s_960s),
  eval_surv_cohort("GSE96058_TNBC_Ext_OS", g960$clin$os_days[g960_e] / 365.25, g960$clin$os_event[g960_e], s_960e),
  eval_surv_cohort("GSE96058_Basal_OS", g960$clin$os_days[g960_b] / 365.25, g960$clin$os_event[g960_b], s_960b)
))
print(surv_summary)
fwrite(surv_summary, file.path(res_dir, "GMSig4_Prognostic_Validation_Summary.csv"))

# 2. Save individual patient-level risk scores + metadata for plotting Figure 3 KM & ROC & Heatmaps
tcga_df <- data.frame(
  sampleID = tcga$clin$sampleID[idx_t],
  OS_yr = tcga$clin$OS.time[idx_t] / 365.25,
  OS = tcga$clin$OS[idx_t],
  GMSig4 = s_tcga,
  RiskGroup = ifelse(s_tcga > median(s_tcga), "High", "Low"),
  Age = as.numeric(tcga$clin$age_at_initial_pathologic_diagnosis[idx_t]),
  Stage_raw = as.character(tcga$clin$pathologic_stage[idx_t]),
  T_raw = as.character(tcga$clin$pathologic_T[idx_t]),
  N_raw = as.character(tcga$clin$pathologic_N[idx_t]),
  Z_tcga,
  stringsAsFactors = FALSE
)
tcga_df$Stage_III_IV <- ifelse(grepl("Stage III|Stage IV", tcga_df$Stage_raw), 1,
                        ifelse(grepl("Stage I|Stage II", tcga_df$Stage_raw), 0, NA))
tcga_df$Node_Pos <- ifelse(grepl("^N[123]", tcga_df$N_raw), 1, ifelse(grepl("^N0", tcga_df$N_raw), 0, NA))
tcga_df$T_large  <- ifelse(grepl("^T[34]", tcga_df$T_raw), 1, ifelse(grepl("^T[12]", tcga_df$T_raw), 0, NA))

# Univariate & Multivariable Cox in TCGA-TNBC
cat("\n=== Multivariable Cox in TCGA-TNBC ===\n")
mv_tcga <- coxph(Surv(OS_yr, OS) ~ GMSig4 + Age + Stage_III_IV + Node_Pos, data = tcga_df)
mv_t_sum <- summary(mv_tcga)
print(mv_t_sum)

# Multivariable Cox + Nomogram in GSE96058_TNBC_Ext (n=242, cleanest full clinical variables!)
g960_df <- data.frame(
  sampleID = g960$clin$sample_title[g960_e],
  OS_yr = g960$clin$os_days[g960_e] / 365.25,
  OS = g960$clin$os_event[g960_e],
  GMSig4 = s_960e,
  RiskGroup = ifelse(s_960e > median(s_960e), "High", "Low"),
  Age = as.numeric(g960$clin$age_at_diagnosis[g960_e]),
  Tumor_mm = as.numeric(g960$clin$tumor_size[g960_e]),
  Node_group = as.character(g960$clin$lymph_node_group[g960_e]),
  NHG = as.character(g960$clin$nhg[g960_e]),
  Chemo = as.character(g960$clin$chemo_treated[g960_e]),
  scale(t(g960$expr[sig_genes, g960_e, drop = FALSE])),
  stringsAsFactors = FALSE
)
g960_df$Node_Pos <- ifelse(g960_df$Node_group %in% c("NodeNegative", "N0"), 0,
                    ifelse(!is.na(g960_df$Node_group) & g960_df$Node_group != "", 1, NA))
g960_df$Chemo_Yes <- suppressWarnings(as.numeric(g960_df$Chemo))

cat("\n=== Multivariable Cox in GSE96058 Extended TNBC (n=", nrow(g960_df), ") ===\n", sep = "")
mv_g960 <- coxph(Surv(OS_yr, OS) ~ GMSig4 + Age + Tumor_mm + Node_Pos, data = g960_df)
mv_9_sum <- summary(mv_g960)
print(mv_9_sum)

# Export Multivariable Cox table for both TCGA and GSE96058
extract_cox_dt <- function(fit, coh_name) {
  s <- summary(fit)
  data.table(
    cohort = coh_name,
    variable = rownames(s$coefficients),
    n = fit$n, events = fit$nevent,
    HR = s$conf.int[, "exp(coef)"],
    LCI = s$conf.int[, "lower .95"],
    UCI = s$conf.int[, "upper .95"],
    pval = s$coefficients[, "Pr(>|z|)"]
  )
}
mv_all <- rbind(extract_cox_dt(mv_tcga, "TCGA_TNBC"), extract_cox_dt(mv_g960, "GSE96058_TNBC_Ext"))
fwrite(mv_all, file.path(res_dir, "GMSig4_Multivariable_Cox_Results.csv"))

# Nomogram C-index comparison: Clinical only vs GMSig4+Clinical (TCGA-TNBC & GSE96058)
tcga_cc <- tcga_df[complete.cases(tcga_df[, c("OS_yr", "OS", "GMSig4", "Age", "Stage_III_IV", "Node_Pos")]), ]
f_clin_t <- coxph(Surv(OS_yr, OS) ~ Age + Stage_III_IV + Node_Pos, data = tcga_cc)
f_nomo_t <- coxph(Surv(OS_yr, OS) ~ GMSig4 + Age + Stage_III_IV + Node_Pos, data = tcga_cc)
g960_cc <- g960_df[complete.cases(g960_df[, c("OS_yr", "OS", "GMSig4", "Age", "Tumor_mm", "Node_Pos")]), ]
f_clin_9 <- coxph(Surv(OS_yr, OS) ~ Age + Tumor_mm + Node_Pos, data = g960_cc)
f_nomo_9 <- coxph(Surv(OS_yr, OS) ~ GMSig4 + Age + Tumor_mm + Node_Pos, data = g960_cc)

nomo_cmp <- data.table(
  cohort = c("TCGA_TNBC", "TCGA_TNBC", "GSE96058_TNBC_Ext", "GSE96058_TNBC_Ext"),
  model = c("Clinical_Only", "Nomogram_GMSig4+Clinical", "Clinical_Only", "Nomogram_GMSig4+Clinical"),
  n = c(nrow(tcga_cc), nrow(tcga_cc), nrow(g960_cc), nrow(g960_cc)),
  C_index = c(summary(f_clin_t)$concordance[1], summary(f_nomo_t)$concordance[1],
              summary(f_clin_9)$concordance[1], summary(f_nomo_9)$concordance[1]),
  LRT_add_GMSig4_p = c(NA_real_, anova(f_clin_t, f_nomo_t)[2, "Pr(>|Chi|)"],
                       NA_real_, anova(f_clin_9, f_nomo_9)[2, "Pr(>|Chi|)"])
)
cat("\n=== Nomogram Incremental Prognostic Value ===\n")
print(nomo_cmp)
fwrite(nomo_cmp, file.path(res_dir, "GMSig4_Nomogram_Comparison.csv"))

# Save per-patient survival data for clean plotting in Figures 3 & 4
g58_df <- data.frame(
  sampleID = rownames(g58$clin)[g58_idx],
  OS_yr = g58$clin$os_days[g58_idx] / 365.25, OS = g58$clin$os_event[g58_idx],
  MFS_yr = g58$clin$mfs_days[g58_idx] / 365.25, MFS = g58$clin$meta_event[g58_idx],
  GMSig4 = s_g58, RiskGroup = ifelse(s_g58 > median(s_g58), "High", "Low"),
  scale(t(g58$expr[sig_genes, g58_idx, drop = FALSE])),
  stringsAsFactors = FALSE
)
fwrite(tcga_df, file.path(res_dir, "GMSig4_PatientScores_TCGA_TNBC.csv"))
fwrite(g58_df,  file.path(res_dir, "GMSig4_PatientScores_GSE58812.csv"))
fwrite(g960_df, file.path(res_dir, "GMSig4_PatientScores_GSE96058_TNBC_Ext.csv"))
saveRDS(list(weights = w_ridge, uni_betas = w_uni, nomo_tcga = f_nomo_t, nomo_g960 = f_nomo_9),
        file.path(res_dir, "GMSig4_Model_Objects.rds"))
cat("SAVED GMSig4 clinical validation & nomogram results!\n")
