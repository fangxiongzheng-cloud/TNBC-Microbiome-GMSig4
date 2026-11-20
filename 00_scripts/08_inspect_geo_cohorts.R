project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
})

data_dir <- file.path(project_root, "01_data", "TCGA_GEO")

# 1. Inspect GSE58812 & GSE103091
for (fn in c("GSE58812_series_matrix.txt.gz", "GSE103091_series_matrix.txt.gz")) {
  fp <- file.path(data_dir, fn)
  if (!file.exists(fp)) next
  cat("\n=== Inspecting", fn, "===\n")
  gset <- getGEO(filename = fp, getGPL = FALSE)
  pd <- pData(gset)
  ex <- exprs(gset)
  cat(sprintf("Samples=%d, Features=%d\n", ncol(ex), nrow(ex)))
  char_cols <- grep("characteristics_ch1", colnames(pd), value = TRUE)
  cat("Characteristic cols:\n")
  for (cc in char_cols) {
    cat(" ", cc, "->", paste(head(unique(pd[[cc]]), 4), collapse = " | "), "\n")
  }
  cat("Feature annotation cols:", paste(colnames(fData(gset)), collapse = ", "), "\n")
}

# 2. Inspect GSE194040 (I-SPY2-990) & GSE173839 (I-SPY2 Durvalumab)
for (fn in c("GSE194040-GPL30493_series_matrix.txt.gz", "GSE194040-GPL20078_series_matrix.txt.gz", "GSE173839_series_matrix.txt.gz")) {
  fp <- file.path(data_dir, fn)
  if (!file.exists(fp)) next
  cat("\n=== Inspecting", fn, "===\n")
  gset <- getGEO(filename = fp, getGPL = FALSE)
  pd <- pData(gset)
  cat(sprintf("Samples=%d\n", nrow(pd)))
  char_cols <- grep("characteristics_ch1|title|source_name", colnames(pd), value = TRUE)
  for (cc in char_cols) {
    cat(" ", cc, "->", paste(head(unique(pd[[cc]]), 5), collapse = " | "), "\n")
  }
}

# 3. Inspect GSE194040 geneLevel_n988 & GSE173839 biomarkers
f_bio <- file.path(data_dir, "GSE173839_biomarkers.csv.gz")
if (file.exists(f_bio) && file.info(f_bio)$size > 100) {
  bio <- fread(f_bio, data.table = FALSE)
  cat("\n=== GSE173839 biomarkers ===\n")
  cat("Rows:", nrow(bio), "Cols:", paste(colnames(bio), collapse = ", "), "\n")
  print(head(bio, 3))
}
