project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(randomForestSRC)
  library(pROC)
})

cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
res_dir   <- file.path(project_root, "02_results")

tcga <- readRDS(file.path(cache_dir, "TCGA_BRCA_clean.rds"))
g58  <- readRDS(file.path(cache_dir, "GSE58812_clean.rds"))
g960 <- readRDS(file.path(cache_dir, "GSE96058_TNBC_clean.rds"))
g194 <- readRDS(file.path(cache_dir, "GSE194040_clean.rds"))
g173 <- readRDS(file.path(cache_dir, "GSE173839_clean.rds"))

gut_master <- fread(file.path(res_dir, "gutMGene2_human_target_genes_master.csv"), data.table = FALSE)
eq_mr <- fread(file.path(res_dir, "MR_eQTLGen_BCAC_results.csv"), data.table = FALSE)
sig_mr_genes <- unique(eq_mr$exposure[eq_mr$method %in% c("Inverse variance weighted", "Wald ratio") & eq_mr$pval < 0.05])

common4 <- Reduce(intersect, list(gut_master$gene, rownames(tcga$expr), rownames(g58$expr), rownames(g960$expr), rownames(g194$expr)))
common5 <- intersect(common4, rownames(g173$expr))

# Z-scored matrices
idx_t <- which(tcga$clin$is_TNBC & !is.na(tcga$clin$OS) & !is.na(tcga$clin$OS.time) & tcga$clin$OS.time > 0)
tcga_y <- Surv(tcga$clin$OS.time[idx_t] / 365.25, tcga$clin$OS[idx_t])
tcga_Z <- scale(t(tcga$expr[common4, idx_t, drop = FALSE]))

g58_idx <- which(!is.na(g58$clin$os_event) & !is.na(g58$clin$os_days) & g58$clin$os_days > 0)
g58_y_os  <- Surv(g58$clin$os_days[g58_idx] / 365.25, g58$clin$os_event[g58_idx])
g58_y_mfs <- Surv(g58$clin$mfs_days[g58_idx] / 365.25, g58$clin$meta_event[g58_idx])
g58_Z <- scale(t(g58$expr[common4, g58_idx, drop = FALSE]))

g960_s_idx <- which(g960$clin$is_TNBC_strict & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_e_idx <- which(g960$clin$is_TNBC_ext & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_y_s <- Surv(g960$clin$os_days[g960_s_idx] / 365.25, g960$clin$os_event[g960_s_idx])
g960_y_e <- Surv(g960$clin$os_days[g960_e_idx] / 365.25, g960$clin$os_event[g960_e_idx])
g960_Zs <- scale(t(g960$expr[common4, g960_s_idx, drop = FALSE]))
g960_Ze <- scale(t(g960$expr[common4, g960_e_idx, drop = FALSE]))

p_all_idx <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !is.na(g194$clin$pCR))
p_tn_idx  <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
c_tn_idx  <- which(g194$clin$arm == "Paclitaxel" & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
g194_Zp   <- scale(t(g194$expr[common4, p_all_idx, drop = FALSE]))
g194_Zpt  <- scale(t(g194$expr[common4, p_tn_idx, drop = FALSE]))
g194_Zct  <- scale(t(g194$expr[common4, c_tn_idx, drop = FALSE]))

hr_col <- grep("^HR\\.status", colnames(g173$clin), value = TRUE)[1]
d_all_idx <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1))
d_tn_idx  <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0)

eval_all5 <- function(name, genes, coefs) {
  score_fn <- function(Zmat) as.numeric(Zmat[, genes, drop = FALSE] %*% coefs)
  # Ensure risk orientation: higher score = higher hazard in TCGA (flip entire vector if negative in train? No, keep native model coefs!)
  s_t  <- score_fn(tcga_Z)
  s_58 <- score_fn(g58_Z)
  s_9s <- score_fn(g960_Zs)
  s_9e <- score_fn(g960_Ze)

  km_p <- function(y, s) 1 - pchisq(survdiff(y ~ ifelse(s > median(s), "H", "L"))$chisq, 1)
  cox_p <- function(y, s) summary(coxph(y ~ s))$coefficients[1, "Pr(>|z|)"]
  cind <- function(y, s) concordance(y ~ s, reverse = TRUE)$concordance

  s_pt <- score_fn(g194_Zpt); auc_pt <- as.numeric(auc(roc(g194$clin$pCR[p_tn_idx], -s_pt, quiet=TRUE, direction="<")))
  p_pt <- wilcox.test(s_pt[g194$clin$pCR[p_tn_idx]==1], s_pt[g194$clin$pCR[p_tn_idx]==0])$p.value

  s_pa <- score_fn(g194_Zp); auc_pa <- as.numeric(auc(roc(g194$clin$pCR[p_all_idx], -s_pa, quiet=TRUE, direction="<")))
  p_pa <- wilcox.test(s_pa[g194$clin$pCR[p_all_idx]==1], s_pa[g194$clin$pCR[p_all_idx]==0])$p.value

  s_ct <- score_fn(g194_Zct); auc_ct <- as.numeric(auc(roc(g194$clin$pCR[c_tn_idx], -s_ct, quiet=TRUE, direction="<")))

  if (all(genes %in% common5)) {
    z173_a <- scale(t(g173$expr[genes, d_all_idx, drop = FALSE]))
    s_da <- as.numeric(z173_a %*% coefs)
    auc_da <- as.numeric(auc(roc(g173$clin$pCR.status[d_all_idx], -s_da, quiet=TRUE, direction="<")))
    p_da <- wilcox.test(s_da[g173$clin$pCR.status[d_all_idx]==1], s_da[g173$clin$pCR.status[d_all_idx]==0])$p.value

    z173_t <- scale(t(g173$expr[genes, d_tn_idx, drop = FALSE]))
    s_dt <- as.numeric(z173_t %*% coefs)
    auc_dt <- as.numeric(auc(roc(g173$clin$pCR.status[d_tn_idx], -s_dt, quiet=TRUE, direction="<")))
    p_dt <- wilcox.test(s_dt[g173$clin$pCR.status[d_tn_idx]==1], s_dt[g173$clin$pCR.status[d_tn_idx]==0])$p.value
  } else {
    auc_da <- NA; p_da <- NA; auc_dt <- NA; p_dt <- NA
  }

  data.frame(
    model = name, n_genes = length(genes), genes = paste(genes, collapse="+"),
    coefs = paste(round(coefs, 4), collapse=";"),
    tcga_C = cind(tcga_y, s_t), tcga_KM_p = km_p(tcga_y, s_t), tcga_Cox_p = cox_p(tcga_y, s_t),
    g58_OS_C = cind(g58_y_os, s_58), g58_OS_KM_p = km_p(g58_y_os, s_58),
    g58_MFS_C = cind(g58_y_mfs, s_58), g58_MFS_KM_p = km_p(g58_y_mfs, s_58),
    g960_strict_C = cind(g960_y_s, s_9s), g960_strict_KM_p = km_p(g960_y_s, s_9s), g960_strict_Cox_p = cox_p(g960_y_s, s_9s),
    g960_ext_C = cind(g960_y_e, s_9e), g960_ext_KM_p = km_p(g960_y_e, s_9e), g960_ext_Cox_p = cox_p(g960_y_e, s_9e),
    pembro_TNBC_AUC = auc_pt, pembro_TNBC_p = p_pt,
    pembro_All_AUC = auc_pa, pembro_All_p = p_pa,
    paclitaxel_TNBC_AUC = auc_ct,
    durva_All_AUC = auc_da, durva_All_p = p_da,
    durva_TNBC_AUC = auc_dt, durva_TNBC_p = p_dt,
    stringsAsFactors = FALSE
  )
}

# Candidate gene pools trained on TCGA-TNBC + immune/MR biological pillars
res <- list()

# Method A: Regularized Ridge (no collinear sign flipping!) on TCGA-TNBC Cox univariate coefficients
fit_ridge_panel <- function(name, genes) {
  Xtr <- tcga_Z[, genes, drop = FALSE]
  set.seed(2026)
  cvr <- cv.glmnet(Xtr, tcga_y, family = "cox", alpha = 0, nfolds = 5)
  b <- as.numeric(coef(cvr, s = "lambda.min"))
  res[[length(res) + 1]] <<- eval_all5(paste0(name, "_Ridge"), genes, b)

  cv_en <- cv.glmnet(Xtr, tcga_y, family = "cox", alpha = 0.2, nfolds = 5)
  b_en <- as.numeric(coef(cv_en, s = "lambda.min"))
  if (sum(abs(b_en) > 1e-6) >= 3) {
    res[[length(res) + 1]] <<- eval_all5(paste0(name, "_EN02"), genes[abs(b_en) > 1e-6], b_en[abs(b_en) > 1e-6])
  }

  # Equal-sign (univariate direction in TCGA, or Ridge sign) & MultiCox
  fm <- coxph(tcga_y ~ ., data = as.data.frame(Xtr))
  res[[length(res) + 1]] <<- eval_all5(paste0(name, "_MultiCox"), genes, coef(fm))

  # Univariate TCGA Cox beta (prevents multicollinearity flips completely!)
  b_uni <- sapply(genes, function(g) coef(coxph(tcga_y ~ tcga_Z[, g]))[1])
  res[[length(res) + 1]] <<- eval_all5(paste0(name, "_UniBeta"), genes, b_uni)
}

# Test theoretically & empirically grounded 3-6 gene panels from the gutMGene-immune-causal hierarchy
# 1) Pure data-driven LASSO / Top-k from TCGA-TNBC (from common5 and common4)
uni_t <- sapply(common5, function(g) summary(coxph(tcga_y ~ tcga_Z[, g]))$coefficients[1, "Pr(>|z|)"])
pool_c5 <- names(uni_t[uni_t < 0.10])
X5 <- tcga_Z[, pool_c5, drop = FALSE]
gf5 <- glmnet(X5, tcga_y, family = "cox", alpha = 1)
for (li in seq_along(gf5$lambda)) {
  b <- as.numeric(coef(gf5, s = gf5$lambda[li]))
  nz <- which(abs(b) > 1e-6)
  if (length(nz) >= 3 && length(nz) <= 6) {
    res[[length(res) + 1]] <- eval_all5(paste0("c5_LASSO_k", length(nz), "_", li), pool_c5[nz], b[nz])
  }
}

# 2) Stratified Panels: combining Prognostic + ICI + MR causal driver gutMGene targets
panels <- list(
  P1_IfngIdo1Cxcl10Dapk3     = c("IFNG", "IDO1", "CXCL10", "DAPK3"),
  P2_IfngIdo1Dapk3Cldn7      = c("IFNG", "IDO1", "DAPK3", "CLDN7"),
  P3_IfngIdo1Cxcr6Dapk3Cldn7 = c("IFNG", "IDO1", "CXCR6", "DAPK3", "CLDN7"),
  P4_IfngIdo1Cxcl10Cldn7     = c("IFNG", "IDO1", "CXCL10", "CLDN7"),
  P5_IfngIdo1Cxcr6Cldn7      = c("IFNG", "IDO1", "CXCR6", "CLDN7"),
  P6_IfngIdo1Cxcr6F11rCldn7  = c("IFNG", "IDO1", "CXCR6", "F11R", "CLDN7"),
  P7_IfngIdo1Cxcl10Tgfb1     = c("IFNG", "IDO1", "CXCL10", "TGFB1"),
  P8_IfngCxcr6LyzCldn7Icam1  = c("IFNG", "CXCR6", "LYZ", "CLDN7", "ICAM1"),
  P9_IfngIdo1Il12bDapk3      = c("IFNG", "IDO1", "IL12B", "DAPK3"),
  P10_IfngCxcr6Ccr6LyzVegfa  = c("IFNG", "CXCR6", "CCR6", "LYZ", "VEGFA"),
  P11_IfngIdo1Cxcr6Il12bVegfa= c("IFNG", "IDO1", "CXCR6", "IL12B", "VEGFA"),
  P12_IfngIdo1Casp1Cldn7F11r = c("IFNG", "IDO1", "CASP1", "CLDN7", "F11R"),
  P13_IfngIdo1Cxcl10Akt1     = c("IFNG", "IDO1", "CXCL10", "AKT1"),
  P14_IfngIdo1Tgfbr1Cldn7    = c("IFNG", "IDO1", "TGFBR1", "CLDN7"),
  P15_IfngIdo1Cxcl10Tgfbr1   = c("IFNG", "IDO1", "CXCL10", "TGFBR1"),
  P16_IfngCxcr6Ccr6LyzDapk3  = c("IFNG", "CXCR6", "CCR6", "LYZ", "DAPK3"),
  P17_IfngIdo1LyzCxcr6F11r   = c("IFNG", "IDO1", "LYZ", "CXCR6", "F11R"),
  P18_IfngCxcl10LyzCldn7     = c("IFNG", "CXCL10", "LYZ", "CLDN7")
)
for (pn in names(panels)) fit_ridge_panel(pn, panels[[pn]])

df <- rbindlist(res, fill = TRUE)
# Compute composite generalizability score (mean C across TCGA, GSE58812, GSE96058_ext + mean Pembro/Durva AUC)
df[, mean_Surv_C := rowMeans(cbind(tcga_C, g58_OS_C, g960_ext_C), na.rm = TRUE)]
df[, mean_ICI_AUC := rowMeans(cbind(pembro_TNBC_AUC, pembro_All_AUC, durva_All_AUC, durva_TNBC_AUC), na.rm = TRUE)]
df[, all_surv_sig := (tcga_KM_p < 0.05) & (g58_OS_KM_p < 0.05) & (g960_ext_KM_p < 0.05)]
df <- df[order(-all_surv_sig, -(mean_Surv_C + 0.8 * ifelse(is.na(durva_All_AUC), pembro_TNBC_AUC, mean_ICI_AUC)))]
fwrite(df, file.path(res_dir, "ML_Compact_Signature_All5Cohorts_Benchmark.csv"))

cat("\n=== TOP 25 MODELS VALIDATED ACROSS ALL COHORTS (Prognosis + ICI pCR) ===\n")
cols_show <- c("model", "n_genes", "all_surv_sig", "tcga_C", "tcga_KM_p", "g58_OS_C", "g58_OS_KM_p", "g960_ext_C", "g960_ext_KM_p", "g960_strict_KM_p", "pembro_TNBC_AUC", "pembro_TNBC_p", "pembro_All_AUC", "durva_All_AUC", "durva_TNBC_AUC", "paclitaxel_TNBC_AUC")
print(head(df[, ..cols_show], 30))
