project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(GEOquery)
})

data_dir  <- file.path(project_root, "01_data", "TCGA_GEO")
cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
dir.create(cache_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== 1. Parse GSE96058 Clinical Metadata (GPL11154 + GPL18573) ===\n")
parse_gse96058_pd <- function(fp, plat) {
  lines <- readLines(gzfile(fp), warn = FALSE)
  sub_lines <- lines[grepl("^!Sample_", lines)]
  # Parse title, geo_accession, and characteristics_ch1 lines
  get_row <- function(prefix) {
    hits <- sub_lines[grepl(paste0("^", prefix, "\t"), sub_lines)]
    lapply(hits, function(h) {
      flds <- strsplit(h, "\t")[[1]][-1]
      gsub("^\"|\"$", "", flds)
    })
  }
  titles <- get_row("!Sample_title")[[1]]
  gsms   <- get_row("!Sample_geo_accession")[[1]]
  chars  <- get_row("!Sample_characteristics_ch1")
  df <- data.frame(gsm = gsms, sample_title = titles, platform = plat, stringsAsFactors = FALSE)
  for (ch in chars) {
    kv <- strsplit(ch, ": ", fixed = TRUE)
    keys <- unique(sapply(kv, function(x) if (length(x) >= 1) trimws(x[1]) else ""))
    keys <- keys[keys != ""]
    if (length(keys) == 1) {
      k <- gsub("[^a-zA-Z0-9]+", "_", tolower(keys[1]))
      vals <- sapply(kv, function(x) if (length(x) >= 2) trimws(paste(x[-1], collapse=": ")) else NA_character_)
      df[[k]] <- vals
    }
  }
  df
}

pd1 <- parse_gse96058_pd(file.path(data_dir, "GSE96058-GPL11154_series_matrix.txt.gz"), "GPL11154")
pd2 <- parse_gse96058_pd(file.path(data_dir, "GSE96058-GPL18573_series_matrix.txt.gz"), "GPL18573")
all_cols <- union(colnames(pd1), colnames(pd2))
for (c in setdiff(all_cols, colnames(pd1))) pd1[[c]] <- NA
for (c in setdiff(all_cols, colnames(pd2))) pd2[[c]] <- NA
pd <- rbind(pd1[, all_cols], pd2[, all_cols])
cat("Combined GSE96058 clinical samples:", nrow(pd), "\n")
cat("Column names:\n")
print(colnames(pd))

# Filter primary tumors (exclude repl samples if any duplicate titles; repl titles end with 'repl')
pd$is_repl <- grepl("repl", pd$sample_title, ignore.case = TRUE)
cat("Replicates:", sum(pd$is_repl), "| Non-replicates:", sum(!pd$is_repl), "\n")

pd$er_num   <- suppressWarnings(as.numeric(pd$er_status))
pd$pgr_num  <- suppressWarnings(as.numeric(pd$pgr_status))
pd$her2_num <- suppressWarnings(as.numeric(pd$her2_status))
pd$er_sgc   <- suppressWarnings(as.numeric(pd$er_prediction_sgc))
pd$pgr_sgc  <- suppressWarnings(as.numeric(pd$pgr_prediction_sgc))
pd$her2_sgc <- suppressWarnings(as.numeric(pd$her2_prediction_sgc))

# Definition 1: Strict clinical IHC TNBC (ER=0 & PgR=0 & HER2=0)
pd$is_TNBC_strict <- !pd$is_repl & !is.na(pd$er_num) & pd$er_num == 0 &
                     !is.na(pd$pgr_num) & pd$pgr_num == 0 & !is.na(pd$her2_num) & pd$her2_num == 0
# Definition 2: Extended TNBC (clinical 0, or if clinical NA then sgc prediction == 0)
er_comb   <- ifelse(!is.na(pd$er_num), pd$er_num, pd$er_sgc)
pgr_comb  <- ifelse(!is.na(pd$pgr_num), pd$pgr_num, pd$pgr_sgc)
her2_comb <- ifelse(!is.na(pd$her2_num), pd$her2_num, pd$her2_sgc)
pd$is_TNBC_ext <- !pd$is_repl & !is.na(er_comb) & er_comb == 0 &
                  !is.na(pgr_comb) & pgr_comb == 0 & !is.na(her2_comb) & her2_comb == 0

pd$os_days  <- suppressWarnings(as.numeric(pd$overall_survival_days))
pd$os_event <- suppressWarnings(as.numeric(pd$overall_survival_event))

cat("Strict clinical TNBC:", sum(pd$is_TNBC_strict),
    "(valid OS:", sum(pd$is_TNBC_strict & !is.na(pd$os_days) & pd$os_days > 0 & !is.na(pd$os_event)),
    "events:", sum(pd$os_event[pd$is_TNBC_strict & !is.na(pd$os_days) & pd$os_days > 0], na.rm=TRUE), ")\n")
cat("Extended TNBC (incl SGC imputed NAs):", sum(pd$is_TNBC_ext),
    "(valid OS:", sum(pd$is_TNBC_ext & !is.na(pd$os_days) & pd$os_days > 0 & !is.na(pd$os_event)),
    "events:", sum(pd$os_event[pd$is_TNBC_ext & !is.na(pd$os_days) & pd$os_days > 0], na.rm=TRUE), ")\n")
cat("PAM50 Basal:", sum(!pd$is_repl & pd$pam50_subtype == "Basal", na.rm = TRUE), "\n")

cat("\n=== 2. Load GSE96058 Transformed Expression Matrix (564MB gz -> subset TNBC & Basal + all genes) ===\n")
# Keep TNBC_ext + Basal samples so matrix is fast & lightweight while preserving all TNBC + Basal controls!
keep_pd <- pd[pd$is_TNBC_ext | (!pd$is_repl & pd$pam50_subtype == "Basal"), ]
ex_dt <- fread(file.path(data_dir, "GSE96058_gene_expression_n3273.csv.gz"), data.table = FALSE)
cat("Raw GSE96058 matrix dims:", dim(ex_dt), "\n")
genes <- ex_dt[[1]]
# Match column names to sample_title
common_t <- intersect(keep_pd$sample_title, colnames(ex_dt)[-1])
cat("Matched TNBC/Basal samples in expr:", length(common_t), "/", nrow(keep_pd), "\n")

keep_pd <- keep_pd[match(common_t, keep_pd$sample_title), ]
mat <- as.matrix(ex_dt[, common_t, drop = FALSE])
rownames(mat) <- genes
rm(ex_dt); gc()

saveRDS(list(expr = mat, clin = keep_pd, all_clin = pd), file.path(cache_dir, "GSE96058_TNBC_clean.rds"))
cat("Saved GSE96058_TNBC_clean.rds:", paste(dim(mat), collapse="x"), "\n")
