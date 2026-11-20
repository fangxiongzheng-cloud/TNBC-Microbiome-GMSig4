project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
})

res_dir <- file.path(project_root, "02_results")
iv_mib  <- fread(file.path(res_dir, "MiBioGen_clumped_IVs_1e5.csv"))
eq_iv   <- fread(file.path(res_dir, "eQTLGen_candidate_instruments.csv"))
eq_rep  <- fread(file.path(res_dir, "MR_eQTLGen_5GWAS_CrossSubtype_Long.csv"))

# Replicated / key causal taxa from 5 GWAS
causal_taxa <- c(
  "genus.Dorea.id.1997",
  "genus.Desulfovibrio.id.3173",
  "genus.LachnospiraceaeNK4A136group.id.11319",
  "genus..Eubacteriumxylanophilumgroup.id.14375",
  "genus.FamilyXIIIUCG001.id.11294",
  "genus.Parabacteroides.id.954",
  "genus.Marvinbryantia.id.2005",
  "genus.Sellimonas.id.14369",
  "genus.Oxalobacter.id.2978",
  "genus.Anaerostipes.id.1991",
  "family.Ruminococcaceae.id.2050",
  "order.Lactobacillales.id.1800"
)

# Replicated causal + signature genes in eQTLGen
target_genes <- c("PEMT", "CCR6", "TCF7L2", "MUC1", "NLRP6", "FOXP3", "MAGI3",
                  "ICAM1", "FFAR2", "ST6GAL1", "LDHA", "TGFBR1", "CXCR6",
                  "CASP1", "DAPK3", "IDO1", "IFNG", "CXCL10", "ENO1", "TGFB1")
eq_map <- unique(eq_iv[exposure %in% target_genes, .(gene = exposure, id.exposure)])
print(eq_map)

sub_iv <- iv_mib[bac %in% causal_taxa]
exp_mib <- data.frame(
  SNP = sub_iv$rsID,
  beta.exposure = as.numeric(sub_iv$beta),
  se.exposure = as.numeric(sub_iv$SE),
  effect_allele.exposure = toupper(sub_iv$alt),
  other_allele.exposure = toupper(sub_iv$ref),
  eaf.exposure = NA_real_,
  pval.exposure = as.numeric(sub_iv$P),
  samplesize.exposure = as.numeric(sub_iv$N),
  exposure = sub_iv$bac,
  id.exposure = sub_iv$bac,
  stringsAsFactors = FALSE
)

med_cache <- file.path(res_dir, "Mediation_Step1_eQTLGen_outcomes_raw.rds")
if (file.exists(med_cache)) {
  out_eq <- readRDS(med_cache)
} else {
  out_eq <- extract_outcome_data(
    snps = unique(exp_mib$SNP),
    outcomes = eq_map$id.exposure,
    proxies = FALSE
  )
  saveRDS(out_eq, med_cache)
}
cat("Extracted Step-1 eQTLGen outcome rows:", nrow(out_eq), "\n")

dat1 <- harmonise_data(exposure_dat = exp_mib, outcome_dat = out_eq, action = 2)
dat1 <- dat1[dat1$mr_keep == TRUE, ]
mr1 <- as.data.table(mr(dat1, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median")))
mr1 <- merge(mr1, eq_map, by.x = "id.outcome", by.y = "id.exposure", all.x = TRUE)
fwrite(mr1, file.path(res_dir, "MR_Mediation_Step1_Taxa_to_Gene_all.csv"))

sig1 <- mr1[method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05][order(pval)]
cat("\n=== Step-1 Significant Causal Links: Gut Taxa -> Host Gene eQTL (p < 0.05) ===\n")
print(sig1[, .(taxon = exposure, gene, nsnp, b1 = b, se1 = se, p1 = pval)])

# Combine Step 1 (Taxon -> Gene) + Step 2 (Gene -> BCAC/ER-/TNBC) + Total Effect (Taxon -> BCAC/ER-/TNBC)
mib_all <- fread(file.path(res_dir, "MR_MiBioGen_5GWAS_CrossSubtype_Long.csv"))
med_rows <- list()
for (i in seq_len(nrow(sig1))) {
  tx <- sig1$exposure[i]
  gn <- sig1$gene[i]
  b1 <- sig1$b[i]; se1 <- sig1$se[i]; p1 <- sig1$pval[i]
  sub2 <- eq_rep[exposure == gn & pval < 0.05]
  if (nrow(sub2) == 0) sub2 <- eq_rep[exposure == gn & pval < 0.10]
  for (j in seq_len(nrow(sub2))) {
    oc <- sub2$outcome_id[j]
    b2 <- log(sub2$OR[j]); se2 <- (log(sub2$OR_uci[j]) - log(sub2$OR_lci[j])) / (2 * 1.96); p2 <- sub2$pval[j]
    tot <- mib_all[exposure == tx & outcome_id == oc]
    b_tot <- if (nrow(tot) > 0) log(tot$OR[1]) else NA_real_
    p_tot <- if (nrow(tot) > 0) tot$pval[1] else NA_real_
    b_med <- b1 * b2
    # Sobel delta-method SE
    se_med <- sqrt(b1^2 * se2^2 + b2^2 * se1^2)
    z_med  <- b_med / se_med
    p_med  <- 2 * pnorm(-abs(z_med))
    med_prop <- ifelse(!is.na(b_tot) && abs(b_tot) > 1e-6, 100 * b_med / b_tot, NA_real_)
    med_rows[[length(med_rows) + 1]] <- data.table(
      taxon = tx, mediator_gene = gn, outcome = oc,
      b_step1 = b1, p_step1 = p1,
      b_step2 = b2, p_step2 = p2,
      b_indirect = b_med, se_indirect = se_med, Sobel_p = p_med,
      b_total = b_tot, p_total = p_tot,
      mediated_pct = med_prop
    )
  }
}
if (length(med_rows) > 0) {
  med_dt <- rbindlist(med_rows)[order(Sobel_p)]
  fwrite(med_dt, file.path(res_dir, "MR_Mediation_2Step_Causal_Pathways.csv"))
  cat("\n=== 2-Step Mediation MR Summary (Gut Taxon -> Gene Expression -> BC/ER-/TNBC) ===\n")
  print(med_dt)
}
