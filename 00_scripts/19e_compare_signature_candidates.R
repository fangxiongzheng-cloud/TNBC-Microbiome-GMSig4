project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(glmnet)
  library(pROC)
})

cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
res_dir   <- file.path(project_root, "02_results")

tcga <- readRDS(file.path(cache_dir, "TCGA_BRCA_clean.rds"))
g58  <- readRDS(file.path(cache_dir, "GSE58812_clean.rds"))
g960 <- readRDS(file.path(cache_dir, "GSE96058_TNBC_clean.rds"))
g194 <- readRDS(file.path(cache_dir, "GSE194040_clean.rds"))
g173 <- readRDS(file.path(cache_dir, "GSE173839_clean.rds"))

# Z-scored cohort sets
idx_t <- which(tcga$clin$is_TNBC & !is.na(tcga$clin$OS) & !is.na(tcga$clin$OS.time) & tcga$clin$OS.time > 0)
tcga_y <- Surv(tcga$clin$OS.time[idx_t] / 365.25, tcga$clin$OS[idx_t])

g58_idx <- which(!is.na(g58$clin$os_event) & !is.na(g58$clin$os_days) & g58$clin$os_days > 0)
g58_y_os  <- Surv(g58$clin$os_days[g58_idx] / 365.25, g58$clin$os_event[g58_idx])
g58_y_mfs <- Surv(g58$clin$mfs_days[g58_idx] / 365.25, g58$clin$meta_event[g58_idx])

g960_s_idx <- which(g960$clin$is_TNBC_strict & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_e_idx <- which(g960$clin$is_TNBC_ext & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_b_idx <- which(g960$clin$pam50_subtype == "Basal" & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_y_s <- Surv(g960$clin$os_days[g960_s_idx] / 365.25, g960$clin$os_event[g960_s_idx])
g960_y_e <- Surv(g960$clin$os_days[g960_e_idx] / 365.25, g960$clin$os_event[g960_e_idx])
g960_y_b <- Surv(g960$clin$os_days[g960_b_idx] / 365.25, g960$clin$os_event[g960_b_idx])

p_all_idx <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !is.na(g194$clin$pCR))
p_tn_idx  <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
c_tn_idx  <- which(g194$clin$arm == "Paclitaxel" & g194$clin$is_TNBC & !is.na(g194$clin$pCR))

hr_col <- grep("^HR\\.status", colnames(g173$clin), value = TRUE)[1]
d_all_idx <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1))
d_tn_idx  <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0)
dc_all_idx <- which(g173$clin$Arm == "control" & g173$clin$pCR.status %in% c(0, 1))

get_z <- function(expr, genes, idx) scale(t(expr[genes, idx, drop = FALSE]))

eval_panel <- function(name, genes, method = "UniBeta") {
  Z_t  <- get_z(tcga$expr, genes, idx_t)
  if (method == "UniBeta") {
    coefs <- sapply(genes, function(g) coef(coxph(tcga_y ~ Z_t[, g]))[1])
  } else if (method == "Ridge") {
    set.seed(2026)
    cvr <- cv.glmnet(Z_t, tcga_y, family = "cox", alpha = 0, nfolds = 5)
    coefs <- as.numeric(coef(cvr, s = "lambda.min"))
  } else {
    coefs <- coef(coxph(tcga_y ~ ., data = as.data.frame(Z_t)))
  }

  sc <- function(expr, idx) as.numeric(get_z(expr, genes, idx) %*% coefs)
  km_p <- function(y, s) 1 - pchisq(survdiff(y ~ ifelse(s > median(s), "H", "L"))$chisq, 1)
  cox_p <- function(y, s) summary(coxph(y ~ s))$coefficients[1, "Pr(>|z|)"]
  cind <- function(y, s) concordance(y ~ s, reverse = TRUE)$concordance
  auc_f <- function(y01, s) as.numeric(auc(roc(y01, -s, quiet = TRUE, direction = "<")))
  wil_f <- function(y01, s) wilcox.test(s[y01 == 1], s[y01 == 0])$p.value

  st  <- sc(tcga$expr, idx_t)
  s58 <- sc(g58$expr, g58_idx)
  s9s <- sc(g960$expr, g960_s_idx)
  s9e <- sc(g960$expr, g960_e_idx)
  s9b <- sc(g960$expr, g960_b_idx)
  spt <- sc(g194$expr, p_tn_idx)
  spa <- sc(g194$expr, p_all_idx)
  sct <- sc(g194$expr, c_tn_idx)

  if (all(genes %in% rownames(g173$expr))) {
    sda <- sc(g173$expr, d_all_idx)
    sdt <- sc(g173$expr, d_tn_idx)
    sdc <- sc(g173$expr, dc_all_idx)
    auc_da <- auc_f(g173$clin$pCR.status[d_all_idx], sda); p_da <- wil_f(g173$clin$pCR.status[d_all_idx], sda)
    auc_dt <- auc_f(g173$clin$pCR.status[d_tn_idx], sdt);  p_dt <- wil_f(g173$clin$pCR.status[d_tn_idx], sdt)
    auc_dc <- auc_f(g173$clin$pCR.status[dc_all_idx], sdc)
  } else {
    auc_da <- NA; p_da <- NA; auc_dt <- NA; p_dt <- NA; auc_dc <- NA
  }

  data.frame(
    model = paste0(name, "_", method), n_genes = length(genes), genes = paste(genes, collapse = "+"),
    coefs = paste(round(coefs, 4), collapse = ";"),
    tcga_C = cind(tcga_y, st), tcga_KM_p = km_p(tcga_y, st), tcga_Cox_p = cox_p(tcga_y, st),
    g58_OS_C = cind(g58_y_os, s58), g58_OS_KM_p = km_p(g58_y_os, s58),
    g58_MFS_C = cind(g58_y_mfs, s58), g58_MFS_KM_p = km_p(g58_y_mfs, s58),
    g960s_C = cind(g960_y_s, s9s), g960s_KM_p = km_p(g960_y_s, s9s),
    g960e_C = cind(g960_y_e, s9e), g960e_KM_p = km_p(g960_y_e, s9e),
    g960b_C = cind(g960_y_b, s9b), g960b_KM_p = km_p(g960_y_b, s9b),
    pembro_TN_AUC = auc_f(g194$clin$pCR[p_tn_idx], spt), pembro_TN_p = wil_f(g194$clin$pCR[p_tn_idx], spt),
    pembro_All_AUC = auc_f(g194$clin$pCR[p_all_idx], spa), pembro_All_p = wil_f(g194$clin$pCR[p_all_idx], spa),
    paclitaxel_TN_AUC = auc_f(g194$clin$pCR[c_tn_idx], sct),
    durva_All_AUC = auc_da, durva_All_p = p_da,
    durva_TN_AUC = auc_dt, durva_TN_p = p_dt, durva_Ctrl_AUC = auc_dc,
    stringsAsFactors = FALSE
  )
}

# Candidate stratified panels combining causal MR genes + prognostic + ICI drivers
cands <- list(
  C1_IfngIdo1Cxcl10Tgfb1     = c("IFNG", "IDO1", "CXCL10", "TGFB1"),
  C2_IfngIdo1Cxcl10Tgfbr1    = c("IFNG", "IDO1", "CXCL10", "TGFBR1"),
  C3_IfngIdo1Cxcl10Tgfb1F11r = c("IFNG", "IDO1", "CXCL10", "TGFB1", "F11R"),
  C4_IfngIdo1Cxcl10Tgfb1Cldn7= c("IFNG", "IDO1", "CXCL10", "TGFB1", "CLDN7"),
  C5_IfngIdo1Cxcl10F11r      = c("IFNG", "IDO1", "CXCL10", "F11R"),
  C6_IfngIdo1Cxcl10Magi3     = c("IFNG", "IDO1", "CXCL10", "MAGI3"),
  C7_IfngIdo1Cxcl10F11rMagi3 = c("IFNG", "IDO1", "CXCL10", "F11R", "MAGI3"),
  C8_IfngIdo1Cxcl10Tgfb1Magi3= c("IFNG", "IDO1", "CXCL10", "TGFB1", "MAGI3"),
  C9_IfngIdo1Cxcr6Tgfb1F11r  = c("IFNG", "IDO1", "CXCR6", "TGFB1", "F11R"),
  C10_IfngIdo1Cxcl10VegfaTgfb1=c("IFNG", "IDO1", "CXCL10", "VEGFA", "TGFB1"),
  C11_IfngIdo1Cxcl10Rac1F11r = c("IFNG", "IDO1", "CXCL10", "RAC1", "F11R"),
  C12_IfngIdo1Cxcr6Magi3Cldn7= c("IFNG", "IDO1", "CXCR6", "MAGI3", "CLDN7"),
  C13_IfngIdo1Cxcl10LyzMagi3 = c("IFNG", "IDO1", "CXCL10", "LYZ", "MAGI3"),
  C14_IfngIdo1F11rMagi3      = c("IFNG", "IDO1", "F11R", "MAGI3"),
  C15_IfngIdo1Tgfb1LyzMagi3  = c("IFNG", "IDO1", "TGFB1", "LYZ", "MAGI3"),
  C16_IfngIdo1Cxcl10Eno1     = c("IFNG", "IDO1", "CXCL10", "ENO1"),
  C17_IfngIdo1Cxcl10Foxp3Magi3=c("IFNG", "IDO1", "CXCL10", "FOXP3", "MAGI3"),
  C18_IfngIdo1Cxcl10Eno1Tgfb1 = c("IFNG", "IDO1", "CXCL10", "ENO1", "TGFB1"),
  C19_IfngIdo1Cxcl10Eno1Tgfbr1= c("IFNG", "IDO1", "CXCL10", "ENO1", "TGFBR1"),
  C20_IfngIdo1Cxcl10Eno1Magi3 = c("IFNG", "IDO1", "CXCL10", "ENO1", "MAGI3"),
  C21_IfngIdo1Cxcl10Eno1Tcf7l2= c("IFNG", "IDO1", "CXCL10", "ENO1", "TCF7L2"),
  C22_IfngIdo1Cxcl10Tgfb1Tcf7l2=c("IFNG", "IDO1", "CXCL10", "TGFB1", "TCF7L2"),
  C23_IfngIdo1Cxcl10Eno1Dapk3 = c("IFNG", "IDO1", "CXCL10", "ENO1", "DAPK3"),
  C24_IfngIdo1Cxcl10Tgfb1Dapk3= c("IFNG", "IDO1", "CXCL10", "TGFB1", "DAPK3"),
  C25_IfngIdo1Cxcl10Eno1Pemt  = c("IFNG", "IDO1", "CXCL10", "ENO1", "PEMT")
)

out <- list()
for (nm in names(cands)) {
  for (m in c("UniBeta", "Ridge")) {
    out[[length(out) + 1]] <- eval_panel(nm, cands[[nm]], m)
  }
}
df <- rbindlist(out)
df[, all_surv_sig := (tcga_KM_p < 0.05) & (g58_OS_KM_p < 0.05) & (g960e_KM_p < 0.05) & (g960b_KM_p < 0.05)]
df[, all5_ici := !is.na(durva_All_AUC)]
fwrite(df[order(-all_surv_sig, -all5_ici, -(tcga_C + g58_OS_C + g960e_C + pembro_TN_AUC + ifelse(is.na(durva_All_AUC), 0, durva_All_AUC)))],
       file.path(res_dir, "Signature_Candidate_Selection_All5Cohorts.csv"))

print(df[order(-all_surv_sig, -all5_ici, -(pembro_TN_AUC + ifelse(is.na(durva_All_AUC), pembro_All_AUC, durva_All_AUC))),
         .(model, all_surv_sig, tcga_C, tcga_KM_p, g58_OS_C, g58_OS_KM_p, g960e_C, g960e_KM_p, g960s_KM_p, pembro_TN_AUC, pembro_All_AUC, durva_All_AUC, durva_TN_AUC, paclitaxel_TN_AUC)])

# Also print published GSE173839 benchmarks for comparison
cat("\n=== Published GSE173839 Durvalumab+Olaparib biomarker AUCs ===\n")
for (col in c("PD1", "PDL1", "T.cells_sig", "TIS_sig", "STAT1_sig", "CD8.T.cells_sig", "Dendritic.cells_sig")) {
  if (col %in% colnames(g173$clin)) {
    v_all <- as.numeric(g173$clin[[col]][d_all_idx])
    v_tn  <- as.numeric(g173$clin[[col]][d_tn_idx])
    a1 <- as.numeric(auc(roc(g173$clin$pCR.status[d_all_idx], v_all, quiet = TRUE, direction = "<")))
    a2 <- as.numeric(auc(roc(g173$clin$pCR.status[d_tn_idx], v_tn, quiet = TRUE, direction = "<")))
    cat(sprintf("  %-22s Durva_All AUC=%.3f | Durva_TNBC AUC=%.3f\n", col, a1, a2))
  }
}
