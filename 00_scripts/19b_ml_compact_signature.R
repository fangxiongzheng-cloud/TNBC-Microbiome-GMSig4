project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(randomForestSRC)
  library(pROC)
  library(timeROC)
})

cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
res_dir   <- file.path(project_root, "02_results")

tcga <- readRDS(file.path(cache_dir, "TCGA_BRCA_clean.rds"))
g58  <- readRDS(file.path(cache_dir, "GSE58812_clean.rds"))
g194 <- readRDS(file.path(cache_dir, "GSE194040_clean.rds"))
g173 <- readRDS(file.path(cache_dir, "GSE173839_clean.rds"))

ev <- fread(file.path(res_dir, "Integrated_gutMGene_MultiEvidence_Ranking.csv"), data.table = FALSE)
gut_master <- fread(file.path(res_dir, "gutMGene2_human_target_genes_master.csv"), data.table = FALSE)

# Common gutMGene genes present in TCGA, GSE58812, GSE194040 (and check GSE173839)
common3 <- Reduce(intersect, list(gut_master$gene, rownames(tcga$expr), rownames(g58$expr), rownames(g194$expr)))
common4 <- intersect(common3, rownames(g173$expr))
cat("gutMGene genes in [TCGA, GSE58812, GSE194040]:", length(common3), "| also in GSE173839:", length(common4), "\n")

# Prepare TCGA-TNBC OS data (train)
idx_t <- which(tcga$clin$is_TNBC & !is.na(tcga$clin$OS) & !is.na(tcga$clin$OS.time) & tcga$clin$OS.time > 0)
tcga_clin <- tcga$clin[idx_t, ]
tcga_y <- Surv(tcga_clin$OS.time / 365.25, tcga_clin$OS)
tcga_X_raw <- t(tcga$expr[common3, idx_t, drop = FALSE])
# Standardize gene expression per cohort (Z-score) for cross-platform comparability
tcga_mu <- colMeans(tcga_X_raw, na.rm = TRUE)
tcga_sd <- apply(tcga_X_raw, 2, sd, na.rm = TRUE)
tcga_sd[tcga_sd == 0] <- 1
tcga_Z <- scale(tcga_X_raw)

# Prepare GSE58812 OS & MFS
g58_clin <- g58$clin
g58_os_idx <- which(!is.na(g58_clin$os_event) & !is.na(g58_clin$os_days) & g58_clin$os_days > 0)
g58_y_os <- Surv(g58_clin$os_days[g58_os_idx] / 365.25, g58_clin$os_event[g58_os_idx])
g58_y_mfs <- Surv(g58_clin$mfs_days[g58_os_idx] / 365.25, g58_clin$meta_event[g58_os_idx])
g58_Z <- scale(t(g58$expr[common3, g58_os_idx, drop = FALSE]))

# Prepare GSE194040 Pembrolizumab & Paclitaxel
p_idx <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !is.na(g194$clin$pCR))
p_tn_idx <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
ctrl_tn_idx <- which(g194$clin$arm == "Paclitaxel" & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
g194_Z_p <- scale(t(g194$expr[common3, p_idx, drop = FALSE]))
g194_Z_ptn <- scale(t(g194$expr[common3, p_tn_idx, drop = FALSE]))
g194_Z_ctn <- scale(t(g194$expr[common3, ctrl_tn_idx, drop = FALSE]))

# Prepare GSE173839 Durvalumab/Olaparib
d_idx <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1))
hr_col <- grep("^HR\\.status", colnames(g173$clin), value = TRUE)[1]
d_tn_idx <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0)

# Helper to evaluate any fixed gene-weight signature across all cohorts
eval_sig <- function(name, genes, coefs) {
  # TCGA-TNBC
  s_tcga <- as.numeric(tcga_Z[, genes, drop = FALSE] %*% coefs)
  c_tcga <- concordance(tcga_y ~ s_tcga, reverse = TRUE)$concordance
  grp_t <- ifelse(s_tcga > median(s_tcga), "High", "Low")
  p_km_t <- 1 - pchisq(survdiff(tcga_y ~ grp_t)$chisq, 1)
  fit_t <- coxph(tcga_y ~ s_tcga)
  p_cox_t <- summary(fit_t)$coefficients[1, "Pr(>|z|)"]

  # GSE58812 OS & MFS
  s_58 <- as.numeric(g58_Z[, genes, drop = FALSE] %*% coefs)
  c_58 <- concordance(g58_y_os ~ s_58, reverse = TRUE)$concordance
  grp_58 <- ifelse(s_58 > median(s_58), "High", "Low")
  p_km_58 <- 1 - pchisq(survdiff(g58_y_os ~ grp_58)$chisq, 1)
  p_cox_58 <- summary(coxph(g58_y_os ~ s_58))$coefficients[1, "Pr(>|z|)"]
  c_58_mfs <- concordance(g58_y_mfs ~ s_58, reverse = TRUE)$concordance
  p_km_58_mfs <- 1 - pchisq(survdiff(g58_y_mfs ~ grp_58)$chisq, 1)

  # GSE194040 Pembrolizumab (lower risk score should predict pCR=1 -> direction=">" or AUC of -score)
  s_ptn <- as.numeric(g194_Z_ptn[, genes, drop = FALSE] %*% coefs)
  auc_ptn <- as.numeric(auc(roc(g194$clin$pCR[p_tn_idx], -s_ptn, quiet = TRUE, direction = "<")))
  p_ptn <- wilcox.test(s_ptn[g194$clin$pCR[p_tn_idx] == 1], s_ptn[g194$clin$pCR[p_tn_idx] == 0])$p.value

  s_pall <- as.numeric(g194_Z_p[, genes, drop = FALSE] %*% coefs)
  auc_pall <- as.numeric(auc(roc(g194$clin$pCR[p_idx], -s_pall, quiet = TRUE, direction = "<")))
  p_pall <- wilcox.test(s_pall[g194$clin$pCR[p_idx] == 1], s_pall[g194$clin$pCR[p_idx] == 0])$p.value

  s_ctn <- as.numeric(g194_Z_ctn[, genes, drop = FALSE] %*% coefs)
  auc_ctn <- as.numeric(auc(roc(g194$clin$pCR[ctrl_tn_idx], -s_ctn, quiet = TRUE, direction = "<")))
  p_ctn <- wilcox.test(s_ctn[g194$clin$pCR[ctrl_tn_idx] == 1], s_ctn[g194$clin$pCR[ctrl_tn_idx] == 0])$p.value

  # GSE173839 Durva (only if all genes present)
  if (all(genes %in% rownames(g173$expr))) {
    z173_a <- scale(t(g173$expr[genes, d_idx, drop = FALSE]))
    s_d <- as.numeric(z173_a %*% coefs)
    auc_d <- as.numeric(auc(roc(g173$clin$pCR.status[d_idx], -s_d, quiet = TRUE, direction = "<")))
    p_d <- wilcox.test(s_d[g173$clin$pCR.status[d_idx] == 1], s_d[g173$clin$pCR.status[d_idx] == 0])$p.value

    z173_t <- scale(t(g173$expr[genes, d_tn_idx, drop = FALSE]))
    s_dt <- as.numeric(z173_t %*% coefs)
    auc_dt <- as.numeric(auc(roc(g173$clin$pCR.status[d_tn_idx], -s_dt, quiet = TRUE, direction = "<")))
    p_dt <- wilcox.test(s_dt[g173$clin$pCR.status[d_tn_idx] == 1], s_dt[g173$clin$pCR.status[d_tn_idx] == 0])$p.value
  } else {
    auc_d <- NA; p_d <- NA; auc_dt <- NA; p_dt <- NA
  }

  data.frame(
    model = name, n_genes = length(genes), genes = paste(genes, collapse="+"),
    tcga_C = c_tcga, tcga_KM_p = p_km_t, tcga_Cox_p = p_cox_t,
    g58_OS_C = c_58, g58_OS_KM_p = p_km_58, g58_OS_Cox_p = p_cox_58,
    g58_MFS_C = c_58_mfs, g58_MFS_KM_p = p_km_58_mfs,
    pembro_TNBC_AUC = auc_ptn, pembro_TNBC_p = p_ptn,
    pembro_All_AUC = auc_pall, pembro_All_p = p_pall,
    paclitaxel_TNBC_AUC = auc_ctn, paclitaxel_TNBC_p = p_ctn,
    durva_All_AUC = auc_d, durva_All_p = p_d,
    durva_TNBC_AUC = auc_dt, durva_TNBC_p = p_dt,
    stringsAsFactors = FALSE
  )
}

# Univariate Cox screen in TCGA-TNBC (pure training set)
uni_t <- sapply(common3, function(g) {
  fit <- coxph(tcga_y ~ tcga_Z[, g])
  summary(fit)$coefficients[1, "Pr(>|z|)"]
})
uni_hr <- sapply(common3, function(g) {
  fit <- coxph(tcga_y ~ tcga_Z[, g])
  summary(fit)$coefficients[1, "exp(coef)"]
})

pool_p05 <- names(uni_t[uni_t < 0.05])
pool_p10 <- names(uni_t[uni_t < 0.10])
pool_c4_p10 <- intersect(pool_p10, common4)
pool_imm_mr <- intersect(common3, c("IFNG", "IDO1", "CXCL10", "CXCR6", "LYZ", "ICAM1", "CCR6", "FOXP3", "ST6GAL1", "MAGI3", "LDHA", "TGFBR1", "CTSD", "DAPK3", "CLDN7", "CASP1", "MYD88", "CCL2"))

cat("TCGA-TNBC uni p<0.05 genes (", length(pool_p05), "):", paste(pool_p05, collapse=", "), "\n")
cat("TCGA-TNBC uni p<0.10 & in GSE173839 (", length(pool_c4_p10), "):", paste(pool_c4_p10, collapse=", "), "\n")

results_list <- list()

# Run Lasso / ElasticNet / StepCox / RSF across pools, constraining K <= 6
set.seed(2026)
for (pname in c("p05", "p10", "c4_p10", "imm_mr")) {
  pool <- switch(pname, p05=pool_p05, p10=pool_p10, c4_p10=pool_c4_p10, imm_mr=pool_imm_mr)
  Xtr <- tcga_Z[, pool, drop = FALSE]

  # 1. StepCox (forward/both) capped at 4, 5, 6 genes via RSF importance or univariate p
  ord_g <- names(sort(uni_t[pool]))
  for (k in c(3, 4, 5, 6)) {
    if (length(ord_g) >= k) {
      gk <- ord_g[1:k]
      f <- coxph(tcga_y ~ ., data = as.data.frame(Xtr[, gk, drop = FALSE]))
      results_list[[length(results_list)+1]] <- eval_sig(sprintf("%s_Top%d_MultiCox", pname, k), gk, coef(f))
    }
  }

  # 2. LASSO path -> pick exact lambda giving 4-6 non-zero genes
  cvf <- cv.glmnet(Xtr, tcga_y, family = "cox", alpha = 1, nfolds = 5)
  gf <- glmnet(Xtr, tcga_y, family = "cox", alpha = 1)
  for (li in seq_along(gf$lambda)) {
    b <- as.numeric(coef(gf, s = gf$lambda[li]))
    nz <- which(abs(b) > 1e-6)
    if (length(nz) >= 3 && length(nz) <= 6) {
      key <- sprintf("%s_LASSO_k%d", pname, length(nz))
      if (!any(sapply(results_list, function(r) r$model == key))) {
        results_list[[length(results_list)+1]] <- eval_sig(key, pool[nz], b[nz])
      }
    }
  }

  # 3. RSF variable importance -> Top 4, 5, 6 MultiCox
  rsf <- rfsrc(Surv(OS.time, OS) ~ ., data = data.frame(OS.time = tcga_clin$OS.time, OS = tcga_clin$OS, Xtr), ntree = 500, importance = TRUE)
  vimp <- sort(rsf$importance, decreasing = TRUE)
  for (k in c(4, 5, 6)) {
    if (length(vimp) >= k) {
      gk <- names(vimp)[1:k]
      f <- coxph(tcga_y ~ ., data = as.data.frame(Xtr[, gk, drop = FALSE]))
      results_list[[length(results_list)+1]] <- eval_sig(sprintf("%s_RSF_Top%d_Cox", pname, k), gk, coef(f))
    }
  }

  # 4. Stepwise AIC on Top 10 univariate -> prune if >6
  top10 <- ord_g[1:min(10, length(ord_g))]
  f_full <- coxph(tcga_y ~ ., data = as.data.frame(Xtr[, top10, drop = FALSE]))
  f_step <- suppressWarnings(step(f_full, direction = "both", trace = 0))
  gs <- names(coef(f_step))
  if (length(gs) >= 3 && length(gs) <= 6) {
    results_list[[length(results_list)+1]] <- eval_sig(sprintf("%s_StepAIC_k%d", pname, length(gs)), gs, coef(f_step))
  }
}

df_bench <- rbindlist(results_list, fill = TRUE)
df_bench <- df_bench[order(-g58_OS_C, -pembro_TNBC_AUC)]
fwrite(df_bench, file.path(res_dir, "ML_Compact_Signature_Benchmark.csv"))
cat("\n=== Top 25 Compact (3-6 gene) Models by Cross-Cohort C-index & ICI AUC ===\n")
print(head(df_bench[, .(model, n_genes, tcga_C, tcga_KM_p, g58_OS_C, g58_OS_KM_p, g58_MFS_C, g58_MFS_KM_p, pembro_TNBC_AUC, pembro_TNBC_p, durva_All_AUC, durva_All_p, paclitaxel_TNBC_AUC, genes)], 25))
