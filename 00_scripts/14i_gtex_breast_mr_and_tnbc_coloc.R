project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

suppressPackageStartupMessages({
  library(data.table)
  library(jsonlite)
  library(TwoSampleMR)
  library(ieugwasr)
  library(coloc)
})

data_dir <- file.path(project_root, "01_data")
res_dir  <- file.path(project_root, "02_results")

coords <- fromJSON(file.path(res_dir, "Candidate_Genes_GRCh38_Coords.json"))
ensg_to_sym <- setNames(names(coords), sapply(coords, `[[`, "ensembl_id"))
ensg_to_sym["ENSG00000126262"] <- "FFAR2"
ensg_to_sym["ENSG00000153721"] <- "CNKSR3"

cat("=== 1. Loading GTEx V8 Breast eQTLs (N=396) & TNBC Locus Summary Stats ===\n")
gtex <- fread(file.path(data_dir, "eQTL", "GTEx_V8_Breast_target_genes_500kb.csv"))
gtex[, ensg := sub("\\..*$", "", gene_id)]
gtex[, symbol := ensg_to_sym[ensg]]
gtex <- gtex[!is.na(symbol) & !is.na(rsid) & rsid != "" & !is.na(beta) & !is.na(se) & se > 0 & !is.na(pvalue)]
gtex[, chr_clean := gsub("^chr", "", as.character(chromosome))]
gtex[, eaf := ifelse(!is.na(ac) & !is.na(an) & an > 0, ac / an, maf)]
gtex <- gtex[eaf > 0.01 & eaf < 0.99]
cat("GTEx Breast genes present:", length(unique(gtex$symbol)), "\n")
print(gtex[, .(n_snps = .N, min_p = min(pvalue)), by = symbol][order(min_p)])

tnbc_loc <- fread(file.path(data_dir, "GWAS_MR", "GCST010100_coloc_regions_500kb.csv"))
tnbc_loc[, chr_clean := gsub("^chr", "", as.character(chromosome))]
tnbc_loc <- tnbc_loc[!is.na(beta) & !is.na(standard_error) & standard_error > 0 & !is.na(p_value)]

# 2. Colocalization (coloc.abf) between GTEx Breast eQTL and TNBC GWAS (GCST010100)
cat("\n=== 2. Running coloc.abf: GTEx V8 Breast cis-eQTL x TNBC (GCST010100) ===\n")
coloc_rows <- list()
for (sym in sort(unique(gtex$symbol))) {
  sub_eq <- gtex[symbol == sym]
  # deduplicate by chr+pos keeping lowest p
  setorder(sub_eq, pvalue)
  sub_eq <- sub_eq[!duplicated(paste(chr_clean, position))]
  ch <- sub_eq$chr_clean[1]
  sub_gw <- tnbc_loc[chr_clean == ch]
  setorder(sub_gw, p_value)
  sub_gw <- sub_gw[!duplicated(base_pair_location)]

  m <- merge(sub_eq, sub_gw, by.x = c("chr_clean", "position"), by.y = c("chr_clean", "base_pair_location"))
  if (nrow(m) < 20) next
  # Match alleles cleanly (or keep MAF/p/beta+varbeta aligned)
  m[, snp_key := paste0(chr_clean, ":", position)]
  # coloc does not depend on allele flip when testing variance explained if MAF & varbeta given,
  # though standardizing p + MAF or aligned beta+varbeta is cleanest:
  d1 <- list(
    beta = m$beta.x,
    varbeta = (m$se)^2,
    MAF = pmin(m$eaf, 1 - m$eaf),
    N = 396,
    type = "quant",
    snp = m$snp_key
  )
  d2 <- list(
    beta = m$beta.y,
    varbeta = (m$standard_error)^2,
    MAF = pmin(m$effect_allele_frequency, 1 - m$effect_allele_frequency),
    N = 118987,
    s = 0.15,
    type = "cc",
    snp = m$snp_key
  )
  # Filter valid MAF
  ok <- !is.na(d1$MAF) & d1$MAF > 0.005 & !is.na(d2$MAF) & d2$MAF > 0.005
  if (sum(ok) < 20) next
  d1 <- lapply(d1, function(v) if (length(v) == nrow(m)) v[ok] else v)
  d2 <- lapply(d2, function(v) if (length(v) == nrow(m)) v[ok] else v)

  res_c <- tryCatch(coloc.abf(dataset1 = d1, dataset2 = d2), error = function(e) NULL)
  if (!is.null(res_c)) {
    sm <- res_c$summary
    top_r <- res_c$results[which.max(res_c$results$SNP.PP.H4), ]
    top_rs <- m[snp_key == top_r$snp]$rsid.x[1]
    coloc_rows[[length(coloc_rows) + 1]] <- data.table(
      gene = sym, n_shared_snps = as.integer(sm["nsnps"]),
      min_gtex_p = min(sub_eq$pvalue), min_tnbc_p = min(m$p_value),
      PP.H0 = sm["PP.H0.abf"], PP.H1 = sm["PP.H1.abf"],
      PP.H2 = sm["PP.H2.abf"], PP.H3 = sm["PP.H3.abf"],
      PP.H4 = sm["PP.H4.abf"],
      cond_H4 = sm["PP.H4.abf"] / (sm["PP.H3.abf"] + sm["PP.H4.abf"]),
      lead_snp = top_rs, lead_SNP_PP_H4 = top_r$SNP.PP.H4
    )
  }
}
coloc_dt <- rbindlist(coloc_rows)[order(-PP.H4)]
fwrite(coloc_dt, file.path(res_dir, "Coloc_GTExBreast_TNBC_results.csv"))
print(coloc_dt[, .(gene, n_shared_snps, min_gtex_p, min_tnbc_p, PP.H3, PP.H4, cond_H4, lead_snp)])

# 3. GTEx V8 Breast Tissue cis-eQTL MR -> TNBC (GCST010100 & GCST90454344) + Overall BC (GCST010098)
cat("\n=== 3. Running GTEx Breast cis-eQTL MR -> TNBC ===\n")
iv_list <- list()
for (sym in sort(unique(gtex$symbol))) {
  sub <- gtex[symbol == sym & pvalue < 1e-3]
  if (nrow(sub) == 0) next
  sub[, F_stat := (beta / se)^2]
  sub <- sub[F_stat > 10]
  if (nrow(sub) == 0) next
  # Clump via OpenGWAS EUR LD (fallback to min-p lead IV if clump returns 0 or fails)
  cl <- tryCatch({
    df_c <- data.frame(rsid = sub$rsid, pval = sub$pvalue, stringsAsFactors = FALSE)
    ld_clump(df_c, clump_kb = 500, clump_r2 = 0.05, clump_p = 1e-3, pop = "EUR")
  }, error = function(e) NULL)
  if (!is.null(cl) && nrow(cl) > 0) {
    kept <- sub[rsid %in% cl$rsid]
  } else {
    kept <- sub[which.min(pvalue)]
  }
  iv_list[[sym]] <- kept
}
gtex_ivs <- rbindlist(iv_list)
cat("Clumped GTEx Breast cis-eQTL IVs:", nrow(gtex_ivs), "across", length(unique(gtex_ivs$symbol)), "genes\n")
fwrite(gtex_ivs, file.path(res_dir, "GTEx_Breast_clumped_IVs.csv"))

exp_gx <- data.frame(
  SNP = gtex_ivs$rsid,
  chr.exposure = gtex_ivs$chr_clean,
  pos.exposure = as.numeric(gtex_ivs$position),
  effect_allele.exposure = toupper(gtex_ivs$alt),
  other_allele.exposure = toupper(gtex_ivs$ref),
  eaf.exposure = as.numeric(gtex_ivs$eaf),
  beta.exposure = as.numeric(gtex_ivs$beta),
  se.exposure = as.numeric(gtex_ivs$se),
  pval.exposure = as.numeric(gtex_ivs$pvalue),
  samplesize.exposure = 396,
  exposure = gtex_ivs$symbol,
  id.exposure = gtex_ivs$symbol,
  stringsAsFactors = FALSE
)

# Match IVs against dense GCST010100_coloc_regions_500kb + PureTNBC/OverallBC API or regions if present
out_tn <- data.frame(
  SNP = tnbc_loc$rsid,
  chr.outcome = tnbc_loc$chr_clean,
  pos.outcome = as.numeric(tnbc_loc$base_pair_location),
  effect_allele.outcome = toupper(tnbc_loc$effect_allele),
  other_allele.outcome = toupper(tnbc_loc$other_allele),
  eaf.outcome = as.numeric(tnbc_loc$effect_allele_frequency),
  beta.outcome = as.numeric(tnbc_loc$beta),
  se.outcome = as.numeric(tnbc_loc$standard_error),
  pval.outcome = as.numeric(tnbc_loc$p_value),
  samplesize.outcome = 118987,
  outcome = "TNBC (GCST010100)",
  id.outcome = "GCST010100",
  mr_keep.outcome = TRUE,
  stringsAsFactors = FALSE
)
out_tn <- out_tn[!duplicated(out_tn$SNP) & out_tn$SNP != "", ]

dat_gx <- harmonise_data(exposure_dat = exp_gx, outcome_dat = out_tn, action = 2)
dat_gx <- dat_gx[dat_gx$mr_keep == TRUE, ]
mr_gx  <- as.data.table(mr(dat_gx, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median")))
mr_gx[, `:=`(OR = exp(b), OR_lci = exp(b - 1.96 * se), OR_uci = exp(b + 1.96 * se))]
mr_gx[, fdr := p.adjust(pval, method = "BH"), by = method]
fwrite(mr_gx, file.path(res_dir, "MR_GTExBreast_GCST010100_TNBC_results.csv"))

cat("\n=== Significant GTEx V8 Breast cis-eQTL MR -> TNBC (GCST010100) p < 0.10 ===\n")
print(mr_gx[method %in% c("Inverse variance weighted", "Wald ratio")][order(pval), .(exposure, method, nsnp, OR, OR_lci, OR_uci, pval, fdr)])
