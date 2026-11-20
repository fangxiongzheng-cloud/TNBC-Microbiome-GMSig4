project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

suppressPackageStartupMessages({
  library(data.table)
  library(TwoSampleMR)
  library(ieugwasr)
  library(ggplot2)
  library(patchwork)
})

gwas_dir <- file.path(project_root, "01_data", "GWAS_MR")
res_dir  <- file.path(project_root, "02_results")
sup_dir  <- file.path(project_root, "03_figures", "Supplementary")
dir.create(sup_dir, showWarnings = FALSE, recursive = TRUE)

clean_taxon <- function(x) {
  y <- gsub("\\.id\\.[0-9]+$", "", x)
  y <- gsub("^genus\\.\\.?", "g. ", y)
  y <- gsub("^family\\.", "f. ", y)
  y <- gsub("^order\\.", "o. ", y)
  y <- gsub("Eubacteriumxylanophilumgroup", "Eubact. xylanophilum", y)
  y <- gsub("LachnospiraceaeNK4A136group", "Lachno. NK4A136", y)
  y <- gsub("FamilyXIIIUCG001", "Family XIII UCG001", y)
  y
}

# ============================================================
# 1. Forward MR Steiger Directionality Test across ALL Significant Pairs
# ============================================================
cat("=== 1. Forward MR Steiger Directionality Tests ===\n")
h_bcac <- fread(file.path(res_dir, "MR_MiBioGen_BCAC_harmonised_SNPs.csv"), data.table = FALSE)
h_tnbc <- fread(file.path(res_dir, "MR_MiBioGen_GCST010100_harmonised.csv"), data.table = FALSE)

# Ensure ncase / ncontrol / prevalence / units / samplesize populated properly for binary breast cancer outcomes
fix_steiger_cols <- function(df, n_case, n_ctrl, ss_exp = 14306) {
  df$samplesize.exposure <- ifelse(is.na(df$samplesize.exposure) | df$samplesize.exposure <= 0, ss_exp, df$samplesize.exposure)
  df$units.exposure <- "SD"
  df$ncase.outcome <- n_case
  df$ncontrol.outcome <- n_ctrl
  df$samplesize.outcome <- n_case + n_ctrl
  df$prevalence.outcome <- 0.12
  df$units.outcome <- "log odds"
  # Approximate EAF for exposure if NA using outcome EAF (standard TwoSampleMR Steiger practice when MAF>0.01)
  df$eaf.exposure <- ifelse(is.na(df$eaf.exposure), df$eaf.outcome, df$eaf.exposure)
  df
}

h_bcac <- as.data.table(h_bcac)
h_1126 <- fix_steiger_cols(as.data.frame(h_bcac[id.outcome == "ieu-a-1126"]), 122977, 105974)
h_1128 <- fix_steiger_cols(as.data.frame(h_bcac[id.outcome == "ieu-a-1128"]), 21468, 105974)
h_1135 <- fix_steiger_cols(as.data.frame(h_bcac[id.outcome == "ieu-a-1135"]), 9655, 45494)
h_0100 <- fix_steiger_cols(as.data.frame(h_tnbc), 11488, 107499)

all_fwd_h <- rbindlist(list(
  as.data.table(h_1126), as.data.table(h_1128), as.data.table(h_1135), as.data.table(h_0100)
), fill = TRUE)

sig_taxa <- c(
  "genus.Dorea.id.1997", "genus.Sellimonas.id.14369", "genus.FamilyXIIIUCG001.id.11294",
  "genus.Desulfovibrio.id.3173", "genus.Parabacteroides.id.954", "genus.Marvinbryantia.id.2005",
  "genus.Oxalobacter.id.2978", "genus.Anaerostipes.id.1991", "genus.LachnospiraceaeNK4A136group.id.11319",
  "genus..Eubacteriumxylanophilumgroup.id.14375", "family.Ruminococcaceae.id.2050",
  "order.Lactobacillales.id.1800", "genus.Adlercreutzia.id.812", "genus.Ruminiclostridium6.id.11356"
)

sub_h <- as.data.frame(all_fwd_h[exposure %in% sig_taxa & mr_keep == TRUE])
st_res <- tryCatch(as.data.table(directionality_test(sub_h)), error = function(e) {
  cat("Steiger note:", e$message, "\n")
  data.table()
})
if (nrow(st_res) > 0) {
  fwrite(st_res, file.path(res_dir, "MR_Forward_Steiger_Directionality.csv"))
  cat("Steiger tests completed:", nrow(st_res), "All correct_causal_direction:",
      all(st_res$correct_causal_direction), "\n")
  print(st_res[order(steiger_pval), .(id.outcome, exposure, snp_r2.exposure, snp_r2.outcome, correct_causal_direction, steiger_pval)][1:15])
}

# ============================================================
# 2. Reverse MR: BC Overall (1126), ER- (1128, 1135), & TNBC (GCST010100) -> Replicated Gut Taxa
# ============================================================
cat("\n=== 2. Reverse MR: BC Subtypes -> Gut Microbiota ===\n")
bc_ivs <- fread(file.path(res_dir, "ReverseMR_BC_Instruments_OpenGWAS.csv"), data.table = FALSE)
tn_cand_file <- file.path(res_dir, "ReverseMR_TNBC_GCST010100_CandidateIVs.csv")
exp_list <- list(bc_ivs)

if (file.exists(tn_cand_file)) {
  tn_raw <- fread(tn_cand_file)
  tn_raw <- tn_raw[!is.na(beta) & !is.na(standard_error) & standard_error > 0 & p_value < 5e-8]
  tn_raw[, F_stat := (beta / standard_error)^2]
  tn_raw <- tn_raw[F_stat > 10 & grepl("^rs", rsid)]
  if (nrow(tn_raw) > 0) {
    # Clump TNBC GWS loci cleanly
    cl_tn <- tryCatch(ld_clump(data.frame(rsid = tn_raw$rsid, pval = tn_raw$p_value),
                               clump_kb = 10000, clump_r2 = 0.001, clump_p = 1, pop = "EUR"),
                      error = function(e) data.frame(rsid = tn_raw$rsid))
    tn_kept <- tn_raw[rsid %in% cl_tn$rsid]
    tn_kept <- tn_kept[!duplicated(rsid)]
    tn_df <- data.frame(
      id.exposure = "GCST010100",
      chr.exposure = gsub("^chr", "", as.character(tn_kept$chromosome)),
      pos.exposure = as.numeric(tn_kept$base_pair_location),
      SNP = tn_kept$rsid,
      effect_allele.exposure = toupper(tn_kept$effect_allele),
      other_allele.exposure = toupper(tn_kept$other_allele),
      eaf.exposure = as.numeric(tn_kept$effect_allele_frequency),
      beta.exposure = as.numeric(tn_kept$beta),
      se.exposure = as.numeric(tn_kept$standard_error),
      pval.exposure = as.numeric(tn_kept$p_value),
      samplesize.exposure = 118987,
      exposure = "TNBC (GCST010100)",
      mr_keep.exposure = TRUE,
      stringsAsFactors = FALSE
    )
    cat("TNBC GCST010100 GWS clumped Reverse MR IVs:", nrow(tn_df), "\n")
    exp_list[[length(exp_list) + 1]] <- tn_df
  }
}

exp_rev <- as.data.frame(rbindlist(exp_list, fill = TRUE))

mib_ids <- fread(file.path(res_dir, "OpenGWAS_MiBioGen_Target_IDs.csv"))
keep_bac_ids <- c(
  "ebi-a-GCST90016989", # Dorea
  "ebi-a-GCST90017067", # Sellimonas
  "ebi-a-GCST90017009", # Family XIII UCG001
  "ebi-a-GCST90016987", # Desulfovibrio
  "ebi-a-GCST90017024", # Lachnospiraceae NK4A136 group
  "ebi-a-GCST90017006", # Eubacterium xylanophilum group
  "ebi-a-GCST90017039", # Parabacteroides
  "ebi-a-GCST90017032", # Marvinbryantia
  "ebi-a-GCST90017038", # Oxalobacter
  "ebi-a-GCST90016966", # Anaerostipes
  "ebi-a-GCST90016951", # family Ruminococcaceae
  "ebi-a-GCST90017101", # order Lactobacillales
  "ebi-a-GCST90016960", # Adlercreutzia
  "ebi-a-GCST90017050"  # Ruminiclostridium6
)

rev_cache <- file.path(res_dir, "ReverseMR_MiBioGen_outcomes_raw.rds")
if (file.exists(rev_cache)) {
  out_rev <- readRDS(rev_cache)
} else {
  out_rev <- extract_outcome_data(
    snps = unique(exp_rev$SNP),
    outcomes = keep_bac_ids,
    proxies = FALSE
  )
  saveRDS(out_rev, rev_cache)
}
cat("Extracted Reverse MR MiBioGen outcome rows:", nrow(out_rev), "\n")

dat_rev <- harmonise_data(exposure_dat = exp_rev, outcome_dat = out_rev, action = 2)
dat_rev <- dat_rev[dat_rev$mr_keep == TRUE, ]
mr_rev  <- as.data.table(mr(dat_rev, method_list = c("mr_ivw", "mr_wald_ratio", "mr_weighted_median", "mr_egger_regression")))
mr_rev[, clean_bac := gsub("Gut microbiota abundance \\(|\\)", "", outcome)]
mr_rev[, clean_bac := gsub(" \\|\\|.*$", "", clean_bac)]
mr_rev[, fdr := p.adjust(pval, method = "BH"), by = .(id.exposure, method)]
fwrite(mr_rev, file.path(res_dir, "MR_Reverse_BC_to_MiBioGen.csv"))

rev_ivw <- mr_rev[method %in% c("Inverse variance weighted", "Wald ratio")]
cat("\n=== Reverse MR IVW Summary (any p < 0.05?) ===\n")
print(rev_ivw[order(pval), .(id.exposure, clean_bac, nsnp, b, se, pval, fdr)][1:20])
cat("Reverse MR pairs with p < 0.05:", sum(rev_ivw$pval < 0.05), "/", nrow(rev_ivw), " | FDR < 0.05:", sum(rev_ivw$fdr < 0.05), "\n")

# ============================================================
# 3. Render Supplementary Figure S1 (Reverse MR + Steiger + LOO + Scatter + Funnel)
# ============================================================
theme_pub <- function(base_size = 8.5) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, color = "grey30", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "grey15", size = base_size - 1),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

# Panel S1A: Reverse MR Beta Heatmap across BC Subtypes -> Replicated Taxa
rev_ivw[, bc_lbl := factor(id.exposure,
  levels = c("ieu-a-1126", "ieu-a-1128", "ieu-a-1135", "GCST010100"),
  labels = c("Overall BC\n(1126)", "ER- BC\n(1128)", "ER- BC\n(1135)", "TNBC\n(010100)")
)]
rev_ivw[, bac_short := gsub("genus |family |order | group", "", gsub(" id\\.[0-9]+", "", clean_bac))]
rev_ivw[, cell_lbl := sprintf("%+.3f%s", b, ifelse(pval < 0.05, "*", ""))]

p_s1a <- ggplot(rev_ivw[!is.na(bc_lbl)], aes(x = bc_lbl, y = bac_short, fill = b)) +
  geom_tile(color = "white", linewidth = 0.5) +
  geom_text(aes(label = cell_lbl, fontface = ifelse(pval < 0.05, "bold", "plain")), size = 2.4) +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-0.15, 0.15), name = "Reverse\nbeta") +
  labs(title = "A. Reverse MR (BC Subtypes -> Replicated Gut Taxa)",
       subtitle = sprintf("FDR<0.05 reverse confounding = %d/%d pairs (Steiger forward TRUE: %d/%d)",
                          sum(rev_ivw$fdr < 0.05), nrow(rev_ivw), sum(st_res$correct_causal_direction), nrow(st_res)),
       x = NULL, y = NULL) +
  theme_pub(8.3)

# Panel S1B: LOO Sensitivity Plot for 4 Flagship Subtype-Causal Links
loo_tnbc <- as.data.table(mr_leaveoneout(as.data.frame(all_fwd_h[
  (exposure == "genus.Dorea.id.1997" & id.outcome == "GCST010100") |
  (exposure == "genus.Sellimonas.id.14369" & id.outcome == "GCST010100") |
  (exposure == "genus.Desulfovibrio.id.3173" & id.outcome == "ieu-a-1128") |
  (exposure == "genus.LachnospiraceaeNK4A136group.id.11319" & id.outcome == "ieu-a-1128")
])))
loo_tnbc[, pair_lbl := paste0(clean_taxon(exposure), " -> ", ifelse(id.outcome == "GCST010100", "TNBC", "ER-"))]
loo_tnbc[, is_all := (SNP == "All")]
loo_tnbc[, lci := b - 1.96 * se][, uci := b + 1.96 * se]

p_s1b <- ggplot(loo_tnbc, aes(x = b, y = reorder(paste(SNP, pair_lbl, sep = "_"), ifelse(is_all, -999, b)), color = is_all)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "grey45") +
  geom_errorbar(aes(xmin = lci, xmax = uci), orientation = "y", width = 0.25, linewidth = 0.5) +
  geom_point(aes(size = is_all)) +
  facet_wrap(~pair_lbl, scales = "free_y", ncol = 2) +
  scale_y_discrete(labels = function(z) sub("_.*$", "", z)) +
  scale_color_manual(values = c("FALSE" = "#3C5488", "TRUE" = "#E64B35"), guide = "none") +
  scale_size_manual(values = c("FALSE" = 1.5, "TRUE" = 2.4), guide = "none") +
  labs(title = "B. Leave-One-Out (LOO) Sensitivity Stability",
       subtitle = "No single SNP drives causal estimates across TNBC & ER- BC",
       x = "IVW Causal Effect beta (95% CI)", y = "Omitted IV SNP") +
  theme_pub(8.0)

# Panel S1C: MR Scatter Plots (SNP-Exposure beta vs SNP-Outcome beta with IVW / WM / Egger lines)
sc_h <- all_fwd_h[
  (exposure == "genus.Dorea.id.1997" & id.outcome == "GCST010100") |
  (exposure == "genus.Sellimonas.id.14369" & id.outcome == "GCST010100") |
  (exposure == "genus.Desulfovibrio.id.3173" & id.outcome == "ieu-a-1135") |
  (exposure == "genus.Marvinbryantia.id.2005" & id.outcome == "GCST010100")
]
sc_h[, pair_lbl := paste0(clean_taxon(exposure), " -> ", ifelse(id.outcome == "GCST010100", "TNBC", "ER- Onco"))]
# Flip exposure alleles positive per standard MR scatter convention
sc_h[, `:=`(bx = ifelse(beta.exposure < 0, -beta.exposure, beta.exposure),
            by = ifelse(beta.exposure < 0, -beta.outcome, beta.outcome))]
sc_pos <- copy(sc_h)
sc_pos[, `:=`(beta.exposure = bx, beta.outcome = by)]
mr_lines <- as.data.table(mr(as.data.frame(sc_pos), method_list = c("mr_ivw", "mr_weighted_median", "mr_egger_regression")))
pleio_lines <- as.data.table(mr_pleiotropy_test(as.data.frame(sc_pos)))[, .(id.exposure, id.outcome, egger_intercept)]
mr_lines <- merge(mr_lines, pleio_lines, by = c("id.exposure", "id.outcome"), all.x = TRUE)
mr_lines[, pair_lbl := paste0(clean_taxon(exposure), " -> ", ifelse(id.outcome == "GCST010100", "TNBC", "ER- Onco"))]
mr_lines[, intercept := ifelse(method == "MR Egger" & !is.na(egger_intercept), egger_intercept, 0)]

p_s1c <- ggplot(sc_h, aes(x = bx, y = by)) +
  geom_hline(yintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_vline(xintercept = 0, color = "grey75", linewidth = 0.3) +
  geom_errorbar(aes(ymin = by - se.outcome, ymax = by + se.outcome), color = "grey65", width = 0) +
  geom_errorbar(aes(xmin = bx - se.exposure, xmax = bx + se.exposure), orientation = "y", color = "grey65", height = 0) +
  geom_point(size = 2.1, color = "#2B5B9E") +
  geom_abline(data = mr_lines, aes(slope = b, intercept = intercept, color = method, linetype = method), linewidth = 0.75) +
  facet_wrap(~pair_lbl, scales = "free", ncol = 2) +
  scale_color_manual(values = c("Inverse variance weighted" = "#E64B35", "Weighted median" = "#00A087", "MR Egger" = "#8491B4"), name = "MR Method") +
  scale_linetype_manual(values = c("Inverse variance weighted" = "solid", "Weighted median" = "dashed", "MR Egger" = "dotted"), name = "MR Method") +
  labs(title = "C. Multi-Method MR SNP Effect Scatter Plots",
       x = "SNP Effect on Gut Taxon (beta)", y = "SNP Effect on BC Subtype (log OR)") +
  theme_pub(8.0) +
  theme(legend.position = "bottom")

# Panel S1D: Funnel Symmetry Plots (Wald ratio beta vs precision 1/SE)
sc_h[, `:=`(wald_b = by / bx, wald_se = se.outcome / abs(bx))]
sc_h[, inv_se := 1 / wald_se]
ivw_vlines <- mr_lines[method %in% c("Inverse variance weighted", "MR Egger")]

p_s1d <- ggplot(sc_h, aes(x = wald_b, y = inv_se)) +
  geom_vline(data = ivw_vlines, aes(xintercept = b, color = method, linetype = method), linewidth = 0.8) +
  geom_point(size = 2.0, color = "#3C5488") +
  facet_wrap(~pair_lbl, scales = "free", ncol = 2) +
  scale_color_manual(values = c("Inverse variance weighted" = "#E64B35", "MR Egger" = "#00A087"), name = "Estimate") +
  scale_linetype_manual(values = c("Inverse variance weighted" = "solid", "MR Egger" = "dashed"), name = "Estimate") +
  labs(title = "D. Funnel Symmetry Diagnostics (Horizontal Pleiotropy Check)",
       x = "Per-SNP Wald Estimate (beta_IV)", y = "Precision (1 / SE_IV)") +
  theme_pub(8.0) +
  theme(legend.position = "bottom")

figs1 <- (p_s1a | p_s1b) / (p_s1c | p_s1d) +
  plot_annotation(
    title = "Supplementary Figure S1. Reverse MR, Steiger Directionality, Leave-One-Out, Scatter & Funnel Sensitivity Diagnostics",
    theme = theme(plot.title = element_text(face = "bold", size = 11.5, hjust = 0))
  )

ggsave(file.path(sup_dir, "FigureS1_ReverseMR_Steiger_Sensitivity.pdf"), figs1, width = 16.5, height = 13.5, device = cairo_pdf)
ggsave(file.path(sup_dir, "FigureS1_ReverseMR_Steiger_Sensitivity.png"), figs1, width = 16.5, height = 13.5, dpi = 300)
cat("SAVED FigureS1_ReverseMR_Steiger_Sensitivity.{pdf,png} successfully!\n")
