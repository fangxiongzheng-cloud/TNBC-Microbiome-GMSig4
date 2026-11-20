project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(MRPRESSO)
})

res_dir <- file.path(project_root, "02_results")
mr_dt   <- fread(file.path(res_dir, "MR_MiBioGen_BCAC_all_results.csv"))
het_dt  <- fread(file.path(res_dir, "MR_MiBioGen_BCAC_heterogeneity.csv"))
ple_dt  <- fread(file.path(res_dir, "MR_MiBioGen_BCAC_pleiotropy.csv"))
dat_all <- fread(file.path(res_dir, "MR_MiBioGen_BCAC_harmonised_SNPs.csv"), data.table = FALSE)

sig_ivw <- mr_dt[method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05][order(pval)]
cat("=== All IVW/Wald p < 0.05 taxa across ieu-a-1126/1128/1135:", nrow(sig_ivw), "===\n")

# Merge heterogeneity & Egger intercept & WM & Egger
het_ivw <- het_dt[method == "Inverse variance weighted", .(id.exposure, id.outcome, Q, Q_df, Q_pval)]
ple_sub <- ple_dt[, .(id.exposure, id.outcome, egger_intercept, egger_se = se, egger_pval = pval)]
wm_sub  <- mr_dt[method == "Weighted median", .(id.exposure, id.outcome, WM_OR = OR, WM_pval = pval)]
eg_sub  <- mr_dt[method == "MR Egger", .(id.exposure, id.outcome, Egger_OR = OR, Egger_pval = pval)]

summary_sig <- Reduce(function(x, y) merge(x, y, by = c("id.exposure", "id.outcome"), all.x = TRUE),
                      list(sig_ivw[, .(id.exposure, id.outcome, exposure, outcome, method, nsnp, b, se, OR, OR_lci, OR_uci, pval, fdr)],
                           wm_sub, eg_sub, het_ivw, ple_sub))

# Run MR-PRESSO for each significant pair with nsnp >= 4
presso_list <- vector("list", nrow(summary_sig))
set.seed(42)
for (i in seq_len(nrow(summary_sig))) {
  ex <- summary_sig$id.exposure[i]
  oc <- summary_sig$id.outcome[i]
  sub_dat <- dat_all[dat_all$id.exposure == ex & dat_all$id.outcome == oc, ]
  if (nrow(sub_dat) >= 4) {
    res_p <- tryCatch({
      mr_presso(
        BetaOutcome = "beta.outcome",
        BetaExposure = "beta.exposure",
        SdOutcome = "se.outcome",
        SdExposure = "se.exposure",
        OUTLIERtest = TRUE,
        DISTORTIONtest = TRUE,
        data = sub_dat,
        NbDistribution = 1000,
        SignifThreshold = 0.05
      )
    }, error = function(e) NULL)
    if (!is.null(res_p)) {
      g_p <- res_p$`MR-PRESSO results`$`Global Test`$Pvalue
      main_tbl <- res_p$`Main MR results`
      raw_p <- main_tbl[1, "P-value"]
      cor_b <- if (nrow(main_tbl) >= 2 && !is.na(main_tbl[2, "Causal Estimate"])) main_tbl[2, "Causal Estimate"] else NA_real_
      cor_p <- if (nrow(main_tbl) >= 2 && !is.na(main_tbl[2, "P-value"])) main_tbl[2, "P-value"] else NA_real_
      presso_list[[i]] <- data.table(
        id.exposure = ex, id.outcome = oc,
        PRESSO_global_p = as.character(g_p),
        PRESSO_raw_p = raw_p,
        PRESSO_corrected_OR = ifelse(is.na(cor_b), NA_real_, exp(cor_b)),
        PRESSO_corrected_p = cor_p
      )
    } else {
      presso_list[[i]] <- data.table(id.exposure = ex, id.outcome = oc, PRESSO_global_p = NA_character_, PRESSO_raw_p = NA_real_, PRESSO_corrected_OR = NA_real_, PRESSO_corrected_p = NA_real_)
    }
  } else {
    presso_list[[i]] <- data.table(id.exposure = ex, id.outcome = oc, PRESSO_global_p = NA_character_, PRESSO_raw_p = NA_real_, PRESSO_corrected_OR = NA_real_, PRESSO_corrected_p = NA_real_)
  }
}
presso_dt <- rbindlist(presso_list)
summary_sig <- merge(summary_sig, presso_dt, by = c("id.exposure", "id.outcome"), all.x = TRUE)

# Consistent direction across IVW, WM, Egger
summary_sig[, direction_consistent := (OR > 1 & (is.na(WM_OR) | WM_OR > 1) & (is.na(Egger_OR) | Egger_OR > 1)) |
                                      (OR < 1 & (is.na(WM_OR) | WM_OR < 1) & (is.na(Egger_OR) | Egger_OR < 1))]

fwrite(summary_sig[order(pval)], file.path(res_dir, "MR_MiBioGen_BCAC_significant_sensitivity.csv"))
print(summary_sig[order(pval), .(id.outcome, exposure, nsnp, OR, pval, WM_pval, Q_pval, egger_pval, PRESSO_global_p, direction_consistent)])

# Leave-one-out on all significant pairs
sig_dat <- dat_all[paste(dat_all$id.exposure, dat_all$id.outcome) %in% paste(summary_sig$id.exposure, summary_sig$id.outcome), ]
loo <- mr_leaveoneout(sig_dat)
fwrite(loo, file.path(res_dir, "MR_MiBioGen_BCAC_significant_LOO.csv"))
cat("Saved sensitivity & LOO tables.\n")
