project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

suppressPackageStartupMessages({
  library(data.table)
  library(org.Hs.eg.db)
  library(AnnotationDbi)
  library(TwoSampleMR)
  library(ieugwasr)
  library(coloc)
})

res_dir <- file.path(project_root, "02_results")
master <- fread(file.path(res_dir, "MultiCohort_TNBC_gutMGene_evidence_matrix.csv"), data.table = FALSE)

# Prioritized gutMGene candidate genes:
# Any gene showing:
#   - TCGA-TNBC OS p < 0.10 OR GSE58812 OS p < 0.05
#   - OR Pembrolizumab TNBC p < 0.05 OR Durvalumab TNBC p < 0.05
#   - OR anchored_GMrepo & TNBC_vs_Normal FDR < 0.05 & |cor_CytotoxicImmune_rho| > 0.45
cand <- unique(master$gene[
  (!is.na(master$TNBC_OS_p) & master$TNBC_OS_p < 0.10) |
  (!is.na(master$GSE58812_OS_p) & master$GSE58812_OS_p < 0.05) |
  (!is.na(master$pembro_TNBC_p) & master$pembro_TNBC_p < 0.05) |
  (!is.na(master$durva_TNBC_p) & master$durva_TNBC_p < 0.05) |
  (master$anchored_GMrepo & !is.na(master$FDR_TNBC_vs_Normal) & master$FDR_TNBC_vs_Normal < 0.05 &
     !is.na(master$cor_CytotoxicImmune_rho) & abs(master$cor_CytotoxicImmune_rho) > 0.45)
])
cand <- cand[!is.na(cand) & cand != ""]
cat("Prioritized gutMGene genes for eQTLGen MR + Coloc:", length(cand), "\n")
print(cand)

ens <- mapIds(org.Hs.eg.db, keys = cand, column = "ENSEMBL", keytype = "SYMBOL", multiVals = "first")
ens_df <- data.frame(gene = cand, ensembl = as.character(ens), stringsAsFactors = FALSE)
ens_df <- ens_df[!is.na(ens_df$ensembl), ]
ens_df$eqtl_id <- paste0("eqtl-a-", ens_df$ensembl)

# Check which eqtl-a-ENSG exist in OpenGWAS
ginfo <- ieugwasr::gwasinfo(ens_df$eqtl_id)
cat("Available in OpenGWAS eQTLGen:", nrow(ginfo), "/", nrow(ens_df), "\n")
ens_df <- merge(ens_df, ginfo[, c("id", "nsnp", "sample_size")], by.x = "eqtl_id", by.y = "id")

# Extract eQTL instruments (p < 5e-8; fallback 1e-5 if none, r2 < 0.001, kb = 10000)
exp_eqtl_list <- list()
for (i in seq_len(nrow(ens_df))) {
  eid <- ens_df$eqtl_id[i]
  gname <- ens_df$gene[i]
  ex <- tryCatch(
    suppressMessages(extract_instruments(outcomes = eid, p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)),
    error = function(e) NULL
  )
  if (is.null(ex) || nrow(ex) == 0) {
    ex <- tryCatch(
      suppressMessages(extract_instruments(outcomes = eid, p1 = 1e-5, clump = TRUE, r2 = 0.001, kb = 10000)),
      error = function(e) NULL
    )
  }
  if (!is.null(ex) && nrow(ex) > 0) {
    ex$exposure <- gname
    ex$F_stat <- (ex$beta.exposure / ex$se.exposure)^2
    ex <- ex[ex$F_stat > 10, ]
    if (nrow(ex) > 0) exp_eqtl_list[[gname]] <- ex
  }
}

exp_eqtl <- as.data.frame(rbindlist(exp_eqtl_list, fill = TRUE))
cat("Extracted eQTL instruments:", nrow(exp_eqtl), "across", length(unique(exp_eqtl$exposure)), "genes\n")
fwrite(exp_eqtl, file.path(res_dir, "eQTLGen_candidate_instruments.csv"))

# Extract BCAC outcomes (Overall ieu-a-1126, ER- ieu-a-1128, OncoArray ER- ieu-a-1135)
outcomes <- c("ieu-a-1126", "ieu-a-1128", "ieu-a-1135")
out_eqtl <- extract_outcome_data(snps = unique(exp_eqtl$SNP), outcomes = outcomes, proxies = TRUE, rsq = 0.8)
dat_eqtl <- harmonise_data(exp_eqtl, out_eqtl, action = 2)
dat_eqtl <- dat_eqtl[dat_eqtl$mr_keep == TRUE, ]

mr_eqtl <- as.data.table(mr(dat_eqtl, method_list = c("mr_wald_ratio", "mr_ivw", "mr_weighted_median", "mr_egger_regression")))
mr_eqtl[, OR := exp(b)]
mr_eqtl[, OR_lci := exp(b - 1.96 * se)]
mr_eqtl[, OR_uci := exp(b + 1.96 * se)]
mr_eqtl[, fdr := p.adjust(pval, method = "BH"), by = .(id.outcome, method)]
fwrite(mr_eqtl, file.path(res_dir, "MR_eQTLGen_BCAC_results.csv"))

cat("\n=== Significant eQTLGen -> BCAC MR hits (IVW/Wald p < 0.05) ===\n")
sig_eq <- mr_eqtl[method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05][order(pval)]
print(sig_eq[, .(id.outcome, exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)])

# Run Bayesian Colocalization (coloc.abf) for any gene with MR p < 0.05 using OpenGWAS associations around lead eQTL SNP (+/- 500kb)
sig_genes <- unique(sig_eq$exposure)
if (length(sig_genes) > 0) {
  cat("\n=== Running coloc.abf for", length(sig_genes), "MR-significant eQTL genes ===\n")
  coloc_rows <- list()
  for (g in sig_genes) {
    eid <- ens_df$eqtl_id[ens_df$gene == g][1]
    top_snp <- exp_eqtl[exp_eqtl$exposure == g, ]
    top_snp <- top_snp[which.min(top_snp$pval.exposure), "SNP"][1]
    # Query SNP chr/pos via associations or ieugwasr::variants_rsid
    vinfo <- tryCatch(ieugwasr::variants_rsid(top_snp), error = function(e) NULL)
    if (is.null(vinfo) || nrow(vinfo) == 0) next
    chr <- vinfo$chr[1]; pos <- as.numeric(vinfo$pos[1])
    reg <- sprintf("%s:%d-%d", chr, max(1, pos - 500000), pos + 500000)
    a_eq <- tryCatch(ieugwasr::associations(variants = reg, id = eid), error = function(e) NULL)
    if (is.null(a_eq) || nrow(a_eq) < 20) next
    for (oid in unique(sig_eq[exposure == g]$id.outcome)) {
      a_out <- tryCatch(ieugwasr::associations(variants = reg, id = oid), error = function(e) NULL)
      if (is.null(a_out) || nrow(a_out) < 20) next
      m <- merge(a_eq, a_out, by = "rsid", suffixes = c(".eq", ".out"))
      m <- m[!is.na(m$beta.eq) & !is.na(m$se.eq) & m$se.eq > 0 & !is.na(m$beta.out) & !is.na(m$se.out) & m$se.out > 0, ]
      m <- m[!duplicated(m$rsid), ]
      if (nrow(m) < 20) next
      # Determine case prop for BCAC outcomes
      s_prop <- ifelse(oid == "ieu-a-1126", 122977/228951, ifelse(oid == "ieu-a-1128", 21468/127442, 9655/55149))
      d1 <- list(snp = m$rsid, beta = m$beta.eq, varbeta = m$se.eq^2, N = 31684, type = "quant", sdY = 1)
      d2 <- list(snp = m$rsid, beta = m$beta.out, varbeta = m$se.out^2, type = "cc", s = s_prop)
      res_c <- tryCatch(suppressWarnings(coloc.abf(d1, d2)), error = function(e) NULL)
      if (!is.null(res_c)) {
        sm <- res_c$summary
        coloc_rows[[paste(g, oid)]] <- data.frame(
          gene = g, eqtl_id = eid, outcome_id = oid, region = reg, nsnps = sm["nsnps"],
          PP.H0 = sm["PP.H0.abf"], PP.H1 = sm["PP.H1.abf"], PP.H2 = sm["PP.H2.abf"],
          PP.H3 = sm["PP.H3.abf"], PP.H4 = sm["PP.H4.abf"], stringsAsFactors = FALSE
        )
      }
    }
  }
  if (length(coloc_rows) > 0) {
    df_coloc <- as.data.table(rbindlist(coloc_rows))[order(-PP.H4)]
    fwrite(df_coloc, file.path(res_dir, "Coloc_eQTLGen_BCAC_results.csv"))
    print(df_coloc)
  }
}
