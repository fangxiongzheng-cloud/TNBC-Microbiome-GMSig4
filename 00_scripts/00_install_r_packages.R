# Install the R dependencies into the active user's R package library.
user_libraries <- strsplit(Sys.getenv("R_LIBS_USER"), .Platform$path.sep, fixed = TRUE)[[1L]]
if (length(user_libraries) && nzchar(user_libraries[[1L]])) {
  dir.create(user_libraries[[1L]], recursive = TRUE, showWarnings = FALSE)
  .libPaths(c(user_libraries[[1L]], .libPaths()))
}

options(repos = c(
  CRAN = "https://cloud.r-project.org",
  mrcieu = "https://mrcieu.r-universe.dev"
))

cran <- c(
  "data.table", "ggplot2", "ggrepel", "glmnet", "jsonlite", "patchwork",
  "pROC", "ragg", "randomForestSRC", "rms", "scales", "survival",
  "systemfonts", "timeROC", "vegan", "coloc", "ieugwasr", "TwoSampleMR",
  "msigdbr", "MicrobiomeStat", "BiocManager"
)
missing_cran <- setdiff(cran, rownames(installed.packages()))
if (length(missing_cran)) install.packages(missing_cran)

bioc <- c(
  "AnnotationDbi", "clusterProfiler", "fgsea", "GEOquery", "GSVA", "hgu133plus2.db",
  "limma", "org.Hs.eg.db", "impute", "preprocessCore", "phyloseq", "ANCOMBC"
)
missing_bioc <- setdiff(bioc, rownames(installed.packages()))
if (length(missing_bioc)) BiocManager::install(missing_bioc, version = "3.20", update = FALSE, ask = FALSE)

required <- c(cran, bioc)
remaining <- setdiff(required, rownames(installed.packages()))
if (length(remaining)) stop("Install the missing packages in this R environment: ", paste(remaining, collapse = ", "))
