project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(survival)
  library(pROC)
})

data_dir <- file.path(project_root, "01_data", "TCGA_GEO")
res_dir  <- file.path(project_root, "02_results")

cat("=== Check annotation pkgs ===\n")
cat("hgu133plus2.db:", "hgu133plus2.db" %in% rownames(installed.packages()), "\n")
cat("org.Hs.eg.db:", "org.Hs.eg.db" %in% rownames(installed.packages()), "\n")

cat("\n=== 1. Check GSE194040 (I-SPY2-990) Pembrolizumab & TNBC counts ===\n")
g1 <- pData(getGEO(filename = file.path(data_dir, "GSE194040-GPL30493_series_matrix.txt.gz"), getGPL = FALSE))
g2 <- pData(getGEO(filename = file.path(data_dir, "GSE194040-GPL20078_series_matrix.txt.gz"), getGPL = FALSE))
pd_ispy <- rbind(
  g1[, c("title", "characteristics_ch1", "characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.5", "characteristics_ch1.6")],
  g2[, c("title", "characteristics_ch1", "characteristics_ch1.2", "characteristics_ch1.3", "characteristics_ch1.5", "characteristics_ch1.6")]
)
pd_ispy$patient_id <- sub("^ISPY2_", "", pd_ispy$title)
pd_ispy$HR   <- sub("^hr: ", "", pd_ispy$characteristics_ch1.2)
pd_ispy$HER2 <- sub("^her2: ", "", pd_ispy$characteristics_ch1.3)
pd_ispy$pCR  <- as.numeric(sub("^pcr: ", "", pd_ispy$characteristics_ch1.5))
pd_ispy$arm  <- sub("^arm: ", "", pd_ispy$characteristics_ch1.6)
pd_ispy$is_TNBC <- (pd_ispy$HR == "0" & pd_ispy$HER2 == "0")

cat("Arms across I-SPY2:\n")
print(table(pd_ispy$arm, pd_ispy$is_TNBC))
cat("\nPembrolizumab pCR by TNBC:\n")
pembro <- pd_ispy[grepl("Pembrolizumab", pd_ispy$arm, ignore.case = TRUE), ]
print(table(TNBC = pembro$is_TNBC, pCR = pembro$pCR))

cat("\n=== 2. Screen gutMGene genes for pCR prediction in I-SPY2 Pembrolizumab & Paclitaxel ===\n")
expr_ispy <- fread(file.path(data_dir, "GSE194040_geneLevel_n988.txt.gz"), data.table = FALSE)
cat("GSE194040 geneLevel dims:", dim(expr_ispy), "\n")
cat("First 5 col names:", paste(head(colnames(expr_ispy), 6), collapse=", "), "\n")
cat("First 5 row names:", paste(head(expr_ispy[[1]], 5), collapse=", "), "\n")

genes <- expr_ispy[[1]]
mat <- as.matrix(expr_ispy[, -1])
rownames(mat) <- genes
# Match patient IDs (col names might be patient_id or X100899)
cnames <- sub("^X", "", colnames(mat))
colnames(mat) <- cnames
common_pembro <- intersect(pembro$patient_id, colnames(mat))
cat("Matched Pembrolizumab patients in expr:", length(common_pembro), "\n")

tcga_screen <- fread(file.path(res_dir, "TCGA_TNBC_gutMGene_screen.csv"), data.table = FALSE)
cand <- intersect(tcga_screen$gene, rownames(mat))
cat("gutMGene targets present in GSE194040:", length(cand), "\n")

p_sub <- pembro[match(common_pembro, pembro$patient_id), ]
m_sub <- mat[cand, common_pembro, drop = FALSE]

pcr_res <- lapply(cand, function(g) {
  x <- as.numeric(m_sub[g, ])
  # All Pembrolizumab (HER2- high risk)
  wt_all <- wilcox.test(x[p_sub$pCR == 1], x[p_sub$pCR == 0])
  auc_all <- as.numeric(auc(roc(p_sub$pCR, x, quiet = TRUE)))
  # TNBC only inside Pembrolizumab
  idx_tn <- which(p_sub$is_TNBC)
  wt_tn <- wilcox.test(x[idx_tn][p_sub$pCR[idx_tn] == 1], x[idx_tn][p_sub$pCR[idx_tn] == 0])
  auc_tn <- as.numeric(auc(roc(p_sub$pCR[idx_tn], x[idx_tn], quiet = TRUE)))
  data.frame(
    gene = g,
    pembro_all_mean_pCR1 = mean(x[p_sub$pCR == 1], na.rm=TRUE),
    pembro_all_mean_pCR0 = mean(x[p_sub$pCR == 0], na.rm=TRUE),
    pembro_all_p = wt_all$p.value,
    pembro_all_AUC = auc_all,
    pembro_TNBC_p = wt_tn$p.value,
    pembro_TNBC_AUC = auc_tn,
    stringsAsFactors = FALSE
  )
})
df_pcr <- rbindlist(pcr_res)
df_pcr <- merge(df_pcr, tcga_screen, by = "gene", all.x = TRUE)
df_pcr <- df_pcr[order(pembro_TNBC_p, pembro_all_p)]
fwrite(df_pcr, file.path(res_dir, "ISPY2_Pembrolizumab_gutMGene_pCR_screen.csv"))

cat("\nTop 25 gutMGene genes predicting Pembrolizumab pCR in TNBC:\n")
print(head(df_pcr[, .(gene, anchored_GMrepo, pembro_TNBC_p, pembro_TNBC_AUC, pembro_all_p, pembro_all_AUC, TNBC_OS_HR, TNBC_OS_p, cor_CytotoxicImmune_rho)], 25))
