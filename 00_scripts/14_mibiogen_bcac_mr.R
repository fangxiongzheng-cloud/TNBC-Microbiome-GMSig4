project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(MRPRESSO)
})

gwas_dir <- file.path(project_root, "01_data", "GWAS_MR")
res_dir  <- file.path(project_root, "02_results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== 1. Load MiBioGen TopHits (P < 1e-4) ===\n")
mbg <- fread(file.path(gwas_dir, "MiBioGen_tophits_all.csv"), data.table = FALSE)
cat("Total MiBioGen rows:", nrow(mbg), "Unique taxa:", length(unique(mbg$bac)), "\n")

# Filter to standard locus-wide IV threshold P < 1e-5, valid rsID, beta, SE
mbg_iv <- mbg[!is.na(mbg$P) & mbg$P < 1e-5 & grepl("^rs", mbg$rsID) & !is.na(mbg$beta) & !is.na(mbg$SE) & mbg$SE > 0, ]
mbg_iv$F_stat <- (mbg_iv$beta / mbg_iv$SE)^2
mbg_iv <- mbg_iv[mbg_iv$F_stat > 10, ]
cat("Filtered P < 1e-5 & F > 10 rows:", nrow(mbg_iv), "Taxa:", length(unique(mbg_iv$bac)), "\n")

# Note: In MiBioGen, 'alt' is effect allele and 'ref' is other allele (check: beta corresponds to alt)
# Verify allele orientation and clump per taxon via OpenGWAS ld_clump (batched carefully)
clump_cache <- file.path(res_dir, "MiBioGen_clumped_IVs_1e5.csv")
if (file.exists(clump_cache)) {
  iv_clumped <- fread(clump_cache, data.table = FALSE)
  cat("Loaded cached clumped IVs:", nrow(iv_clumped), "\n")
} else {
  cat("Clumping IVs per taxon via OpenGWAS (EUR, r2=0.001, kb=10000)...\n")
  taxa <- sort(unique(mbg_iv$bac))
  clumped_list <- list()
  for (i in seq_along(taxa)) {
    tx <- taxa[i]
    sub <- mbg_iv[mbg_iv$bac == tx, ]
    if (nrow(sub) == 1) {
      clumped_list[[tx]] <- sub
      next
    }
    inp <- data.frame(rsid = sub$rsID, pval = sub$P, id = tx, stringsAsFactors = FALSE)
    # Deduplicate rsID keeping min pval
    inp <- inp[order(inp$pval), ]
    inp <- inp[!duplicated(inp$rsid), ]
    res_cl <- tryCatch({
      suppressMessages(ieugwasr::ld_clump(inp, clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, pop = "EUR"))
    }, error = function(e) {
      Sys.sleep(1.5)
      tryCatch(suppressMessages(ieugwasr::ld_clump(inp, clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, pop = "EUR")), error = function(e2) NULL)
    })
    if (!is.null(res_cl) && nrow(res_cl) > 0) {
      clumped_list[[tx]] <- sub[sub$rsID %in% res_cl$rsid, ]
      clumped_list[[tx]] <- clumped_list[[tx]][!duplicated(clumped_list[[tx]]$rsID), ]
    }
    if (i %% 25 == 0) cat(sprintf("  Clumped %d / %d taxa...\n", i, length(taxa)))
  }
  iv_clumped <- as.data.frame(rbindlist(clumped_list))
  fwrite(iv_clumped, clump_cache)
  cat("Saved clumped IVs:", nrow(iv_clumped), "across", length(unique(iv_clumped$bac)), "taxa\n")
}

# Format as TwoSampleMR exposure_dat
exp_dat <- data.frame(
  SNP = iv_clumped$rsID,
  beta.exposure = as.numeric(iv_clumped$beta),
  se.exposure = as.numeric(iv_clumped$SE),
  effect_allele.exposure = iv_clumped$alt,
  other_allele.exposure = iv_clumped$ref,
  eaf.exposure = NA_real_,
  pval.exposure = as.numeric(iv_clumped$P),
  samplesize.exposure = as.numeric(iv_clumped$N),
  exposure = iv_clumped$bac,
  id.exposure = iv_clumped$bac,
  F_stat = as.numeric(iv_clumped$F_stat),
  stringsAsFactors = FALSE
)

# 2. Extract outcomes from OpenGWAS: ieu-a-1126 (Overall BC), ieu-a-1128 (BCAC ER-), ieu-a-1135 (OncoArray ER-)
outcomes <- c("ieu-a-1126", "ieu-a-1128", "ieu-a-1135")
all_snps <- unique(exp_dat$SNP)
cat("\n=== 2. Extracting outcomes for", length(all_snps), "unique clumped SNPs ===\n")

out_cache <- file.path(res_dir, "BCAC_outcomes_OpenGWAS_raw.rds")
if (file.exists(out_cache)) {
  out_dat <- readRDS(out_cache)
} else {
  out_dat <- extract_outcome_data(snps = all_snps, outcomes = outcomes, proxies = TRUE, rsq = 0.8)
  saveRDS(out_dat, out_cache)
}
cat("Extracted outcome rows:", nrow(out_dat), "\n")

# 3. Harmonise and run TwoSampleMR + sensitivity tests
dat <- harmonise_data(exposure_dat = exp_dat, outcome_dat = out_dat, action = 2)
dat <- dat[dat$mr_keep == TRUE, ]
fwrite(dat, file.path(res_dir, "MR_MiBioGen_BCAC_harmonised_SNPs.csv"))

cat("Running MR across all (exposure x outcome) pairs...\n")
mr_res <- mr(dat, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression"))
het_res <- mr_heterogeneity(dat)
pleio_res <- mr_pleiotropy_test(dat)

# Add FDR per outcome & method
mr_dt <- as.data.table(mr_res)
mr_dt[, fdr := p.adjust(pval, method = "BH"), by = .(id.outcome, method)]
mr_dt[, OR := exp(b)]
mr_dt[, OR_lci := exp(b - 1.96 * se)]
mr_dt[, OR_uci := exp(b + 1.96 * se)]

fwrite(mr_dt, file.path(res_dir, "MR_MiBioGen_BCAC_all_results.csv"))
fwrite(het_res, file.path(res_dir, "MR_MiBioGen_BCAC_heterogeneity.csv"))
fwrite(pleio_res, file.path(res_dir, "MR_MiBioGen_BCAC_pleiotropy.csv"))

cat("\n=== Significant IVW/Wald results for ER- (ieu-a-1128 & ieu-a-1135) p < 0.05 ===\n")
sig_er <- mr_dt[id.outcome %in% c("ieu-a-1128", "ieu-a-1135") &
                method %in% c("Inverse variance weighted", "Wald ratio") &
                pval < 0.05][order(pval)]
print(sig_er[, .(id.outcome, exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)])

cat("\n=== Significant IVW/Wald results for Overall BC (ieu-a-1126) p < 0.05 ===\n")
sig_all <- mr_dt[id.outcome == "ieu-a-1126" &
                 method %in% c("Inverse variance weighted", "Wald ratio") &
                 pval < 0.05][order(pval)]
print(head(sig_all[, .(id.outcome, exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)], 25))
