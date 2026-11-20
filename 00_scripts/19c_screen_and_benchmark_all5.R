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
cat("gutMGene present in [TCGA, GSE58812, GSE96058, GSE194040]:", length(common4), "| + GSE173839:", length(common5), "\n")
cat("MR-significant eQTLGen genes in common5:", paste(intersect(sig_mr_genes, common5), collapse=", "), "\n")
cat("MR-significant genes missing in GSE173839:", paste(setdiff(intersect(sig_mr_genes, common4), common5), collapse=", "), "\n")

# Prepare cohorts (per-cohort Z-score for transferability)
idx_t <- which(tcga$clin$is_TNBC & !is.na(tcga$clin$OS) & !is.na(tcga$clin$OS.time) & tcga$clin$OS.time > 0)
tcga_y <- Surv(tcga$clin$OS.time[idx_t] / 365.25, tcga$clin$OS[idx_t])
tcga_Z <- scale(t(tcga$expr[common4, idx_t, drop = FALSE]))

g58_idx <- which(!is.na(g58$clin$os_event) & !is.na(g58$clin$os_days) & g58$clin$os_days > 0)
g58_y_os  <- Surv(g58$clin$os_days[g58_idx] / 365.25, g58$clin$os_event[g58_idx])
g58_y_mfs <- Surv(g58$clin$mfs_days[g58_idx] / 365.25, g58$clin$meta_event[g58_idx])
g58_Z <- scale(t(g58$expr[common4, g58_idx, drop = FALSE]))

g960_s_idx <- which(g960$clin$is_TNBC_strict & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_e_idx <- which(g960$clin$is_TNBC_ext & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))
g960_b_idx <- which(g960$clin$pam50_subtype == "Basal" & !is.na(g960$clin$os_days) & g960$clin$os_days > 0 & !is.na(g960$clin$os_event))

g960_y_s <- Surv(g960$clin$os_days[g960_s_idx] / 365.25, g960$clin$os_event[g960_s_idx])
g960_y_e <- Surv(g960$clin$os_days[g960_e_idx] / 365.25, g960$clin$os_event[g960_e_idx])
g960_y_b <- Surv(g960$clin$os_days[g960_b_idx] / 365.25, g960$clin$os_event[g960_b_idx])

g960_Zs <- scale(t(g960$expr[common4, g960_s_idx, drop = FALSE]))
g960_Ze <- scale(t(g960$expr[common4, g960_e_idx, drop = FALSE]))
g960_Zb <- scale(t(g960$expr[common4, g960_b_idx, drop = FALSE]))

# Immunotherapy trials
p_all_idx <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !is.na(g194$clin$pCR))
p_tn_idx  <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
c_tn_idx  <- which(g194$clin$arm == "Paclitaxel" & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
g194_Zp   <- scale(t(g194$expr[common4, p_all_idx, drop = FALSE]))
g194_Zpt  <- scale(t(g194$expr[common4, p_tn_idx, drop = FALSE]))
g194_Zct  <- scale(t(g194$expr[common4, c_tn_idx, drop = FALSE]))

hr_col <- grep("^HR\\.status", colnames(g173$clin), value = TRUE)[1]
d_all_idx <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1))
d_tn_idx  <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0)

# First check individual genes in GSE96058-TNBC OS too!
cat("\n=== Individual gutMGene prognostic performance in GSE96058 (Strict n=143 / Extended n=242 / Basal n=339) ===\n")
g960_rows <- lapply(common4, function(g) {
  f_s <- summary(coxph(g960_y_s ~ g960_Zs[, g]))$coefficients[1, ]
  f_e <- summary(coxph(g960_y_e ~ g960_Ze[, g]))$coefficients[1, ]
  f_b <- summary(coxph(g960_y_b ~ g960_Zb[, g]))$coefficients[1, ]
  f_t <- summary(coxph(tcga_y ~ tcga_Z[, g]))$coefficients[1, ]
  f_5 <- summary(coxph(g58_y_os ~ g58_Z[, g]))$coefficients[1, ]
  data.frame(
    gene = g,
    in_GSE173839 = g %in% common5,
    in_MR = g %in% sig_mr_genes,
    TCGA_HR = f_t["exp(coef)"], TCGA_p = f_t["Pr(>|z|)"],
    G58_HR = f_5["exp(coef)"], G58_p = f_5["Pr(>|z|)"],
    G960s_HR = f_s["exp(coef)"], G960s_p = f_s["Pr(>|z|)"],
    G960e_HR = f_e["exp(coef)"], G960e_p = f_e["Pr(>|z|)"],
    Basal_HR = f_b["exp(coef)"], Basal_p = f_b["Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})
df960 <- rbindlist(g960_rows)
df960 <- df960[order(G960e_p)]
fwrite(df960, file.path(res_dir, "TCGA_GSE58812_GSE96058_Univariate_Cox_AllGenes.csv"))

# Print genes consistent across >=2 survival cohorts (p<0.10 same direction) or key MR genes
cons3 <- df960[((TCGA_p < 0.15) + (G58_p < 0.10) + (G960e_p < 0.10) >= 2) | gene %in% c("IFNG", "IDO1", "CXCL10", "CXCR6", "LYZ", "ICAM1", "CCR6", "ST6GAL1", "MAGI3", "LDHA", "TGFBR1", "FOXP3", "CLDN7", "CTSD", "DAPK3", "F11R")]
print(cons3[order(G960e_p), .(gene, in_GSE173839, in_MR, TCGA_HR, TCGA_p, G58_HR, G58_p, G960s_HR, G960s_p, G960e_HR, G960e_p, Basal_HR, Basal_p)])
