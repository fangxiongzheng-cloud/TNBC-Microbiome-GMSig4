project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
  library(pROC)
})

data_dir <- file.path(project_root, "01_data", "TCGA_GEO")
res_dir  <- file.path(project_root, "02_results")

cat("=== Inspect GSE173839 (Durvalumab + Olaparib) ===\n")
gset <- getGEO(filename = file.path(data_dir, "GSE173839_series_matrix.txt.gz"), getGPL = FALSE)
pd <- pData(gset)
cat("Samples:", nrow(pd), "Features:", nrow(exprs(gset)), "\n")
for (cc in grep("characteristics_ch1|title", colnames(pd), value = TRUE)) {
  cat(" ", cc, "->", paste(head(unique(pd[[cc]]), 5), collapse = " | "), "\n")
}

bio <- fread(file.path(data_dir, "GSE173839_biomarkers.csv.gz"), data.table = FALSE)
cat("\nBiomarkers columns:", paste(colnames(bio), collapse=", "), "\n")
print(head(bio, 4))
if ("Arm" %in% colnames(bio)) print(table(Arm = bio$Arm, pCR = bio$pCR))
if ("arm" %in% colnames(bio)) print(table(arm = bio$arm, pCR = bio$pCR))
if ("HR" %in% colnames(bio)) print(table(HR = bio$HR, pCR = bio$pCR))
if ("Receptor.Subtype" %in% colnames(bio)) print(table(bio$Receptor.Subtype, bio$pCR))
