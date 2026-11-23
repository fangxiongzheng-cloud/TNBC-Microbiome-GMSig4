.libPaths(c("/root/R/tnbc_microbiome_lib", .libPaths()))

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
})

gwas_dir <- "/mnt/d/Bioinformatics/TNBC_Microbiome/01_data/GWAS_MR"
res_dir  <- "/mnt/d/Bioinformatics/TNBC_Microbiome/02_results"

gwas_raw <- fread(file.path(gwas_dir, "GCST010098_OverallBC_IVs_extracted.csv"))
snp_map  <- fread(file.path(res_dir, "IV_rsids_GRCh38_map.csv"))
snp_map[, chr38_clean := as.character(chr38)]
gwas_raw[, chr_clean := gsub("^chr", "", as.character(chromosome))]

gwas_m <- merge(gwas_raw, snp_map[, .(rsid_iv = rsid, chr38_clean, pos38)],
                by.x = c("chr_clean", "base_pair_location"), by.y = c("chr38_clean", "pos38"),
                all.x = FALSE, allow.cartesian = TRUE)
gwas_m[, exact_rs := (rsid == rsid_iv)]
setorder(gwas_m, rsid_iv, -exact_rs, p_value)
gwas_m <- gwas_m[!duplicated(rsid_iv)]

out_ov <- data.frame(
  SNP = gwas_m$rsid_iv,
  chr.outcome = gwas_m$chr_clean,
  pos.outcome = as.numeric(gwas_m$base_pair_location),
  effect_allele.outcome = toupper(gwas_m$effect_allele),
  other_allele.outcome = toupper(gwas_m$other_allele),
  eaf.outcome = as.numeric(gwas_m$effect_allele_frequency),
  beta.outcome = as.numeric(gwas_m$beta),
  se.outcome = as.numeric(gwas_m$standard_error),
  pval.outcome = as.numeric(gwas_m$p_value),
  outcome = "Overall Breast Cancer (GCST010098, Zhang et al. 2020)",
  id.outcome = "GCST010098",
  mr_keep.outcome = TRUE,
  stringsAsFactors = FALSE
)
cat("GCST010098 unique matched SNPs:", nrow(out_ov), "\n")

iv_mib <- fread(file.path(res_dir, "MiBioGen_clumped_IVs_1e5.csv"))
exp_mib <- data.frame(
  SNP = iv_mib$rsID,
  beta.exposure = as.numeric(iv_mib$beta),
  se.exposure = as.numeric(iv_mib$SE),
  effect_allele.exposure = toupper(iv_mib$alt),
  other_allele.exposure = toupper(iv_mib$ref),
  eaf.exposure = NA_real_,
  pval.exposure = as.numeric(iv_mib$P),
  samplesize.exposure = as.numeric(iv_mib$N),
  exposure = iv_mib$bac,
  id.exposure = iv_mib$bac,
  stringsAsFactors = FALSE
)
dat_mib <- harmonise_data(exposure_dat = exp_mib, outcome_dat = out_ov, action = 2)
dat_mib <- dat_mib[dat_mib$mr_keep == TRUE, ]
mr_mib <- as.data.table(mr(dat_mib, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
mr_mib[, fdr := p.adjust(pval, method = "BH"), by = method]
mr_mib[, OR := exp(b)]
mr_mib[, OR_lci := exp(b - 1.96 * se)]
mr_mib[, OR_uci := exp(b + 1.96 * se)]
fwrite(mr_mib, file.path(res_dir, "MR_MiBioGen_GCST010098_results.csv"))

exp_eq <- fread(file.path(res_dir, "eQTLGen_candidate_instruments.csv"), data.table = FALSE)
dat_eq <- harmonise_data(exposure_dat = exp_eq, outcome_dat = out_ov, action = 2)
dat_eq <- dat_eq[dat_eq$mr_keep == TRUE, ]
mr_eq <- as.data.table(mr(dat_eq, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
mr_eq[method %in% c("Inverse variance weighted", "Wald ratio"), fdr := p.adjust(pval, method = "BH")]
mr_eq[, OR := exp(b)]
mr_eq[, OR_lci := exp(b - 1.96 * se)]
mr_eq[, OR_uci := exp(b + 1.96 * se)]
fwrite(mr_eq, file.path(res_dir, "MR_eQTLGen_GCST010098_results.csv"))

cat("=== Significant MiBioGen -> GCST010098 p < 0.05 ===\n")
print(mr_mib[method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05][order(pval), .(exposure, nsnp, OR, pval)])
cat("=== Significant eQTLGen -> GCST010098 p < 0.05 ===\n")
print(mr_eq[method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05][order(pval), .(exposure, nsnp, OR, pval)])
