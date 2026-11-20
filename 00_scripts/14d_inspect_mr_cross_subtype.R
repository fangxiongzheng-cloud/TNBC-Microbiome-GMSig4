project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

res_dir <- file.path(project_root, "02_results")
tnbc <- fread(file.path(res_dir, "MR_MiBioGen_GCST010100_results.csv"))
het  <- fread(file.path(res_dir, "MR_MiBioGen_GCST010100_heterogeneity.csv"))
ple  <- fread(file.path(res_dir, "MR_MiBioGen_GCST010100_pleiotropy.csv"))
bcac <- fread(file.path(res_dir, "MR_MiBioGen_BCAC_all_results.csv"))

cat("=== 1. MiBioGen -> TNBC (GCST010100) IVW/Wald p < 0.05 ===\n")
sig_t <- tnbc[method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05][order(pval)]
m_t <- merge(sig_t[, .(exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)],
             het[method == "Inverse variance weighted", .(exposure, Q_pval)], by = "exposure", all.x = TRUE)
m_t <- merge(m_t, ple[, .(exposure, egger_pval = pval)], by = "exposure", all.x = TRUE)
m_t <- merge(m_t, tnbc[method == "Weighted median", .(exposure, WM_OR = OR, WM_pval = pval)], by = "exposure", all.x = TRUE)
print(m_t[order(pval)])

cat("\n=== 2. Check ER- causal taxa inside TNBC (GCST010100) ===\n")
er_taxa <- unique(bcac[id.outcome %in% c("ieu-a-1128", "ieu-a-1135") &
                       method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05]$exposure)
print(tnbc[exposure %in% er_taxa & method %in% c("Inverse variance weighted", "Wald ratio"),
           .(exposure, nsnp, OR, OR_lci, OR_uci, pval)][order(pval)])
