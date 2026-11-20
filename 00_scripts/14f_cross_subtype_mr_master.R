project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages(library(data.table))

res_dir <- file.path(project_root, "02_results")
pt   <- fread(file.path(res_dir, "MR_MiBioGen_GCST90454344_results.csv"))
tn   <- fread(file.path(res_dir, "MR_MiBioGen_GCST010100_results.csv"))
bc   <- fread(file.path(res_dir, "MR_MiBioGen_BCAC_all_results.csv"))

cat("=== 1. MiBioGen -> Pure TNBC (GCST90454344) IVW/Wald p < 0.05 ===\n")
print(pt[method %in% c("Inverse variance weighted", "Wald ratio") & pval < 0.05][order(pval), .(exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)])

# Cross-check ALL taxa across BCAC 1126, ER- 1128/1135, TNBC 010100, PureTNBC 90454344
all_m <- rbind(
  bc[method %in% c("Inverse variance weighted", "Wald ratio"), .(outcome_id = id.outcome, exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)],
  tn[method %in% c("Inverse variance weighted", "Wald ratio"), .(outcome_id = "GCST010100", exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)],
  pt[method %in% c("Inverse variance weighted", "Wald ratio"), .(outcome_id = "GCST90454344", exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)]
)
wide_p  <- dcast(all_m, exposure ~ outcome_id, value.var = "pval")
wide_or <- dcast(all_m, exposure ~ outcome_id, value.var = "OR")
sig_multi <- all_m[pval < 0.05, .(n_sig = .N, mean_dir = mean(OR > 1)), by = exposure][n_sig >= 2]
cat("\n=== 2. Taxa Replicated (p < 0.05) Across >=2 Breast Cancer / TNBC Subtype GWAS ===\n")
rep_tbl <- merge(sig_multi, wide_p, by = "exposure")
rep_tbl <- merge(rep_tbl, wide_or, by = "exposure", suffixes = c("_p", "_OR"))
print(rep_tbl[order(-n_sig, GCST90454344_p)])
fwrite(all_m, file.path(res_dir, "MR_MiBioGen_5GWAS_CrossSubtype_Long.csv"))
fwrite(rep_tbl, file.path(res_dir, "MR_MiBioGen_5GWAS_Replicated_Summary.csv"))

# Also combine eQTLGen across 5 GWAS
eq_bc <- fread(file.path(res_dir, "MR_eQTLGen_BCAC_results.csv"))
eq_tn <- fread(file.path(res_dir, "MR_eQTLGen_GCST010100_results.csv"))
eq_pt <- fread(file.path(res_dir, "MR_eQTLGen_GCST90454344_results.csv"))
all_eq <- rbind(
  eq_bc[method %in% c("Inverse variance weighted", "Wald ratio"), .(outcome_id = id.outcome, exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)],
  eq_tn[method %in% c("Inverse variance weighted", "Wald ratio"), .(outcome_id = "GCST010100", exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)],
  eq_pt[method %in% c("Inverse variance weighted", "Wald ratio"), .(outcome_id = "GCST90454344", exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)]
)
eq_sig_m <- all_eq[pval < 0.05, .(n_sig = .N), by = exposure][n_sig >= 2]
eq_wp <- dcast(all_eq, exposure ~ outcome_id, value.var = "pval")
eq_wo <- dcast(all_eq, exposure ~ outcome_id, value.var = "OR")
eq_rep <- merge(merge(eq_sig_m, eq_wp, by = "exposure"), eq_wo, by = "exposure", suffixes = c("_p", "_OR"))
cat("\n=== 3. eQTLGen Genes Replicated (p < 0.05) Across >=2 Breast Cancer / TNBC Subtype GWAS ===\n")
print(eq_rep[order(-n_sig)])
fwrite(all_eq, file.path(res_dir, "MR_eQTLGen_5GWAS_CrossSubtype_Long.csv"))
fwrite(eq_rep, file.path(res_dir, "MR_eQTLGen_5GWAS_Replicated_Summary.csv"))
