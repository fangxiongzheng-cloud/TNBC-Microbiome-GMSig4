project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(hgu133plus2.db)
  library(AnnotationDbi)
})

data_dir <- file.path(project_root, "01_data", "TCGA_GEO")
cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== 1. Prepare TCGA-BRCA TNBC / Non-TNBC / Normal ===\n")
clin <- fread(file.path(data_dir, "TCGA_BRCA_clinicalMatrix.tsv"), data.table = FALSE)
surv <- fread(file.path(data_dir, "TCGA_BRCA_survival.tsv"), data.table = FALSE)
er  <- ifelse(is.na(clin$breast_carcinoma_estrogen_receptor_status) | clin$breast_carcinoma_estrogen_receptor_status == "",
              clin$ER_Status_nature2012, clin$breast_carcinoma_estrogen_receptor_status)
pr  <- ifelse(is.na(clin$breast_carcinoma_progesterone_receptor_status) | clin$breast_carcinoma_progesterone_receptor_status == "",
              clin$PR_Status_nature2012, clin$breast_carcinoma_progesterone_receptor_status)
her_ihc  <- clin$lab_proc_her2_neu_immunohistochemistry_receptor_status
her_fish <- clin$lab_procedure_her2_neu_in_situ_hybrid_outcome_type
her <- ifelse(her_ihc %in% c("Positive", "Negative"), her_ihc,
              ifelse(her_fish %in% c("Positive", "Negative"), her_fish, clin$HER2_Final_Status_nature2012))

clin$is_TNBC    <- (clin$sample_type == "Primary Tumor") & (er == "Negative") & (pr == "Negative") & (her == "Negative")
clin$is_Normal  <- (clin$sample_type == "Solid Tissue Normal")
clin$is_NonTNBC <- (clin$sample_type == "Primary Tumor") & !clin$is_TNBC & (er %in% c("Positive", "Negative"))

expr_dt <- fread(file.path(data_dir, "TCGA_BRCA_HiSeqV2.tsv.gz"), data.table = FALSE)
tcga_mat <- as.matrix(expr_dt[, -1])
rownames(tcga_mat) <- expr_dt[[1]]
rm(expr_dt); gc()

common_s <- intersect(colnames(tcga_mat), clin$sampleID)
tcga_mat <- tcga_mat[, common_s, drop = FALSE]
clin <- clin[match(common_s, clin$sampleID), ]
clin <- merge(clin, surv, by.x = "sampleID", by.y = "sample", all.x = TRUE, sort = FALSE)
clin <- clin[match(common_s, clin$sampleID), ]

saveRDS(list(expr = tcga_mat, clin = clin), file.path(cache_dir, "TCGA_BRCA_clean.rds"))
cat("Saved TCGA_BRCA_clean.rds:", paste(dim(tcga_mat), collapse="x"), "TNBC:", sum(clin$is_TNBC), "\n")

cat("=== 2. Prepare GSE58812 (107 TNBC) ===\n")
g58 <- getGEO(filename = file.path(data_dir, "GSE58812_series_matrix.txt.gz"), getGPL = FALSE)
pd58 <- pData(g58)
raw58 <- exprs(g58)
if (max(raw58, na.rm = TRUE) > 50) raw58 <- log2(pmax(raw58, 0) + 1)
sym58 <- mapIds(hgu133plus2.db, keys = rownames(raw58), column = "SYMBOL", keytype = "PROBEID", multiVals = "first")
keep58 <- !is.na(sym58) & sym58 != ""
dt58 <- as.data.table(raw58[keep58, , drop = FALSE])
dt58[, gene := sym58[keep58]]
dt58_agg <- dt58[, lapply(.SD, mean, na.rm = TRUE), by = gene]
ex58 <- as.matrix(dt58_agg[, -1, with = FALSE])
rownames(ex58) <- dt58_agg$gene
pd58$meta_event <- as.numeric(sub("^meta: ", "", pd58$characteristics_ch1.2))
pd58$mfs_days   <- as.numeric(sub("^mfs \\(days\\): ", "", pd58$characteristics_ch1.3))
pd58$os_event   <- as.numeric(sub("^death: ", "", pd58$characteristics_ch1.4))
pd58$os_days    <- as.numeric(sub("^os \\(days\\): ", "", pd58$characteristics_ch1.5))
saveRDS(list(expr = ex58, clin = pd58), file.path(cache_dir, "GSE58812_clean.rds"))
cat("Saved GSE58812_clean.rds:", paste(dim(ex58), collapse="x"), "\n")

cat("=== 3. Prepare GSE194040 (I-SPY2-990) ===\n")
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

expr_ispy <- fread(file.path(data_dir, "GSE194040_geneLevel_n988.txt.gz"), data.table = FALSE)
m194 <- as.matrix(expr_ispy[, -1])
rownames(m194) <- expr_ispy[[1]]
colnames(m194) <- sub("^X", "", colnames(m194))
common194 <- intersect(pd_ispy$patient_id, colnames(m194))
pd_ispy <- pd_ispy[match(common194, pd_ispy$patient_id), ]
m194 <- m194[, common194, drop = FALSE]
saveRDS(list(expr = m194, clin = pd_ispy), file.path(cache_dir, "GSE194040_clean.rds"))
cat("Saved GSE194040_clean.rds:", paste(dim(m194), collapse="x"), "\n")

cat("=== 4. Prepare GSE173839 (Durvalumab/Olaparib I-SPY2) ===\n")
annot <- fread(file.path(data_dir, "GSE194040_ProbeAnnotation.txt.gz"), data.table = FALSE)
durva_raw <- fread(file.path(data_dir, "GSE173839_durva_probeLevel_n105.txt.gz"), header = FALSE, data.table = FALSE)
hdr_ids <- sub("\\.0+$", "", as.character(durva_raw[1, -1]))
durva <- durva_raw[-1, ]
colnames(durva) <- c("ProbeName", hdr_ids)
bio <- fread(file.path(data_dir, "GSE173839_biomarkers.csv.gz"), data.table = FALSE)
m_idx <- match(durva$ProbeName, annot$ProbeName)
gsym <- annot$GeneName[m_idx]
keep173 <- !is.na(gsym) & gsym != "" & gsym != "---"
dmat <- apply(as.matrix(durva[keep173, -1]), 2, as.numeric)
dt_d <- as.data.table(dmat)
dt_d[, gene := gsym[keep173]]
dt_d_agg <- dt_d[, lapply(.SD, mean, na.rm = TRUE), by = gene]
ex173 <- as.matrix(dt_d_agg[, -1, with = FALSE])
rownames(ex173) <- dt_d_agg$gene
bio$id_str <- as.character(bio$ResearchID)
common173 <- intersect(bio$id_str, colnames(ex173))
bio <- bio[match(common173, bio$id_str), ]
ex173 <- ex173[, common173, drop = FALSE]
saveRDS(list(expr = ex173, clin = bio), file.path(cache_dir, "GSE173839_clean.rds"))
cat("Saved GSE173839_clean.rds:", paste(dim(ex173), collapse="x"), "\n")
