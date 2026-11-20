project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(hgu133plus2.db)
  library(AnnotationDbi)
  library(survival)
  library(pROC)
})

data_dir <- file.path(project_root, "01_data", "TCGA_GEO")
res_dir  <- file.path(project_root, "02_results")

tcga_screen <- fread(file.path(res_dir, "TCGA_TNBC_gutMGene_screen.csv"), data.table = FALSE)
ispy_pembro <- fread(file.path(res_dir, "ISPY2_Pembrolizumab_gutMGene_pCR_screen.csv"), data.table = FALSE)

map_gpl570 <- function(mat) {
  probes <- rownames(mat)
  sym <- mapIds(hgu133plus2.db, keys = probes, column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
  keep <- !is.na(sym) & sym != ""
  mat2 <- mat[keep, , drop = FALSE]
  sym2 <- sym[keep]
  dt <- as.data.table(mat2)
  dt[, gene := sym2]
  dt_agg <- dt[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  out <- as.matrix(dt_agg[, -1, with = FALSE])
  rownames(out) <- dt_agg$gene
  out
}

cat("=== 1. Process GSE58812 (107 TNBC) with log2 check ===\n")
g58 <- getGEO(filename = file.path(data_dir, "GSE58812_series_matrix.txt.gz"), getGPL = FALSE)
pd58 <- pData(g58)
raw58 <- exprs(g58)
cat("GSE58812 raw range:", range(raw58, na.rm=TRUE), "quantile 99%:", quantile(raw58, 0.99, na.rm=TRUE), "\n")
if (max(raw58, na.rm = TRUE) > 50) {
  raw58 <- log2(pmax(raw58, 0) + 1)
  cat("Applied log2(x+1), new range:", range(raw58, na.rm=TRUE), "\n")
}
ex58 <- map_gpl570(raw58)

pd58$meta_event <- as.numeric(sub("^meta: ", "", pd58$characteristics_ch1.2))
pd58$mfs_days   <- as.numeric(sub("^mfs \\(days\\): ", "", pd58$characteristics_ch1.3))
pd58$os_event   <- as.numeric(sub("^death: ", "", pd58$characteristics_ch1.4))
pd58$os_days    <- as.numeric(sub("^os \\(days\\): ", "", pd58$characteristics_ch1.5))

cand <- intersect(tcga_screen$gene, rownames(ex58))
res58 <- lapply(cand, function(g) {
  x <- as.numeric(ex58[g, ])
  fit_os <- coxph(Surv(os_days, os_event) ~ x, data = pd58)
  s_os <- summary(fit_os)
  fit_mfs <- coxph(Surv(mfs_days, meta_event) ~ x, data = pd58)
  s_mfs <- summary(fit_mfs)
  data.frame(
    gene = g,
    GSE58812_OS_HR = s_os$coefficients[1, "exp(coef)"],
    GSE58812_OS_p  = s_os$coefficients[1, "Pr(>|z|)"],
    GSE58812_MFS_HR = s_mfs$coefficients[1, "exp(coef)"],
    GSE58812_MFS_p  = s_mfs$coefficients[1, "Pr(>|z|)"],
    stringsAsFactors = FALSE
  )
})
df58 <- rbindlist(res58)

cat("\n=== 2. Process GSE173839 (Durvalumab + Olaparib) if ProbeAnnotation ready ===\n")
annot_fp <- file.path(data_dir, "GSE194040_ProbeAnnotation.txt.gz")
df173 <- data.frame(gene = character(0), durva_all_AUC = numeric(0), durva_all_p = numeric(0), durva_TNBC_AUC = numeric(0), durva_TNBC_p = numeric(0))
if (file.exists(annot_fp) && file.info(annot_fp)$size > 10000) {
  annot <- fread(annot_fp, data.table = FALSE)
  cat("ProbeAnnotation cols:", paste(colnames(annot), collapse=", "), "\n")
  durva_raw <- fread(file.path(data_dir, "GSE173839_durva_probeLevel_n105.txt.gz"), header = FALSE, data.table = FALSE)
  hdr_ids <- sub("\\.0+$", "", as.character(durva_raw[1, -1]))
  durva <- durva_raw[-1, ]
  colnames(durva) <- c("ProbeName", hdr_ids)
  bio   <- fread(file.path(data_dir, "GSE173839_biomarkers.csv.gz"), data.table = FALSE)
  pcol <- "ProbeName"
  gcol <- "GeneName"

  m_idx <- match(durva$ProbeName, annot[[pcol]])
  gsym <- annot[[gcol]][m_idx]
  keep <- !is.na(gsym) & gsym != "" & gsym != "---"
  dmat <- apply(as.matrix(durva[keep, -1]), 2, as.numeric)
  dt_d <- as.data.table(dmat)
  dt_d[, gene := gsym[keep]]
  dt_d_agg <- dt_d[, lapply(.SD, mean, na.rm = TRUE), by = gene]
  ex173 <- as.matrix(dt_d_agg[, -1, with = FALSE])
  rownames(ex173) <- dt_d_agg$gene
  cat("GSE173839 gene-level dims:", dim(ex173), "\n")

  bio_arm <- bio[bio$Arm == "durvalumab/olaparib" & bio$pCR.status %in% c(0, 1), ]
  bio_arm$id_str <- as.character(bio_arm$ResearchID)
  common_ids <- intersect(bio_arm$id_str, colnames(ex173))
  cat("Matched durvalumab/olaparib samples:", length(common_ids), "\n")
  b_sub <- bio_arm[match(common_ids, bio_arm$id_str), ]
  m_sub <- ex173[, common_ids, drop = FALSE]
  hr_col <- grep("^HR\\.status", colnames(b_sub), value = TRUE)[1]
  is_tnbc_d <- b_sub[[hr_col]] == 0

  cand173 <- intersect(tcga_screen$gene, rownames(m_sub))
  res173 <- lapply(cand173, function(g) {
    x <- as.numeric(m_sub[g, ])
    wt_a <- wilcox.test(x[b_sub$pCR.status == 1], x[b_sub$pCR.status == 0])
    auc_a <- as.numeric(auc(roc(b_sub$pCR.status, x, quiet = TRUE)))
    wt_t <- wilcox.test(x[is_tnbc_d][b_sub$pCR.status[is_tnbc_d] == 1], x[is_tnbc_d][b_sub$pCR.status[is_tnbc_d] == 0])
    auc_t <- as.numeric(auc(roc(b_sub$pCR.status[is_tnbc_d], x[is_tnbc_d], quiet = TRUE)))
    data.frame(gene = g, durva_all_AUC = auc_a, durva_all_p = wt_a$p.value, durva_TNBC_AUC = auc_t, durva_TNBC_p = wt_t$p.value, stringsAsFactors = FALSE)
  })
  df173 <- rbindlist(res173)
}

master <- as.data.table(merge(ispy_pembro, df58, by = "gene", all = TRUE))
if (nrow(df173) > 0) master <- as.data.table(merge(master, df173, by = "gene", all = TRUE))
master <- master[order(GSE58812_OS_p)]
fwrite(master, file.path(res_dir, "MultiCohort_TNBC_gutMGene_evidence_matrix.csv"))

cat("\n=== Genes consistently prognostic (same HR direction, p<0.10 in BOTH TCGA-TNBC & GSE58812) ===\n")
cons_surv <- master[!is.na(TNBC_OS_p) & !is.na(GSE58812_OS_p) &
                    TNBC_OS_p < 0.10 & GSE58812_OS_p < 0.10 &
                    ((TNBC_OS_HR > 1 & GSE58812_OS_HR > 1) | (TNBC_OS_HR < 1 & GSE58812_OS_HR < 1))]
print(cons_surv[, .(gene, anchored_GMrepo, TNBC_OS_HR, TNBC_OS_p, GSE58812_OS_HR, GSE58812_OS_p, pembro_TNBC_AUC, pembro_TNBC_p, durva_TNBC_AUC, durva_TNBC_p)])

cat("\n=== Genes predictive in BOTH Immunotherapy Trials (pembro p<0.05 & durva p<0.10) ===\n")
imm_both <- master[!is.na(pembro_all_p) & !is.na(durva_all_p) & (pembro_TNBC_p < 0.05 | pembro_all_p < 0.05) & (durva_TNBC_p < 0.10 | durva_all_p < 0.10)]
print(imm_both[order(pembro_TNBC_p), .(gene, anchored_GMrepo, pembro_TNBC_AUC, pembro_TNBC_p, pembro_all_AUC, pembro_all_p, durva_TNBC_AUC, durva_TNBC_p, durva_all_AUC, durva_all_p, TNBC_OS_p, GSE58812_OS_p)])
