project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(MicrobiomeStat)
  library(vegan)
})

gm_dir <- file.path(project_root, "01_data", "GMrepo")
res_dir <- file.path(project_root, "02_results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

cohorts <- c("PRJNA658160", "PRJNA383849", "PRJNA872152")
all_diff <- list()
alpha_summary <- list()

for (proj in cohorts) {
  fpath <- file.path(gm_dir, paste0("genus_abundance_", proj, ".csv"))
  if (!file.exists(fpath)) next
  dt <- fread(fpath, data.table = FALSE)
  meta_cols <- c("run_id", "project_id", "sample_id", "disease", "term", "host_age", "sex", "BMI", "country")
  meta <- dt[, meta_cols]
  meta$Group <- ifelse(meta$disease == "D001943", "BC", "Health")
  meta$Group <- factor(meta$Group, levels = c("Health", "BC"))
  rownames(meta) <- meta$run_id

  taxa_cols <- setdiff(colnames(dt), meta_cols)
  mat <- as.matrix(dt[, taxa_cols])
  rownames(mat) <- dt$run_id
  mat[is.na(mat)] <- 0

  # Keep genera present in >= 10% of samples in this cohort
  prev <- colMeans(mat > 0)
  keep <- names(prev)[prev >= 0.10]
  mat_filt <- mat[, keep, drop = FALSE]

  cat(sprintf("\n=== Cohort %s: Total=%d (BC=%d, Health=%d), Genera_kept(>=10%% prev)=%d ===\n",
              proj, nrow(meta), sum(meta$Group == "BC"), sum(meta$Group == "Health"), ncol(mat_filt)))

  # Alpha diversity (Shannon & Richness on relative abundances * 1e4 pseudo-counts)
  shannon <- diversity(mat, index = "shannon")
  richness <- rowSums(mat > 0)
  w_shan <- wilcox.test(shannon ~ meta$Group)
  w_rich <- wilcox.test(richness ~ meta$Group)
  alpha_summary[[proj]] <- data.frame(
    project = proj,
    n_BC = sum(meta$Group == "BC"),
    n_Health = sum(meta$Group == "Health"),
    shannon_BC_median = median(shannon[meta$Group == "BC"]),
    shannon_Health_median = median(shannon[meta$Group == "Health"]),
    shannon_p = w_shan$p.value,
    richness_BC_median = median(richness[meta$Group == "BC"]),
    richness_Health_median = median(richness[meta$Group == "Health"]),
    richness_p = w_rich$p.value
  )

  # Beta diversity (Bray-Curtis PERMANOVA)
  bc_dist <- vegdist(mat_filt, method = "bray")
  perm <- adonis2(bc_dist ~ Group, data = meta, permutations = 999)
  alpha_summary[[proj]]$permanova_R2 <- perm$R2[1]
  alpha_summary[[proj]]$permanova_p <- perm$`Pr(>F)`[1]

  # 1) Wilcoxon test + log2FC on relative abundance
  res_list <- lapply(colnames(mat_filt), function(g) {
    x_bc <- mat_filt[meta$Group == "BC", g]
    x_hc <- mat_filt[meta$Group == "Health", g]
    wt <- wilcox.test(x_bc, x_hc)
    data.frame(
      project = proj,
      genus = g,
      mean_BC = mean(x_bc),
      mean_Health = mean(x_hc),
      prev_BC = mean(x_bc > 0),
      prev_Health = mean(x_hc > 0),
      log2FC = log2((mean(x_bc) + 1e-4) / (mean(x_hc) + 1e-4)),
      wilcox_p = wt$p.value,
      stringsAsFactors = FALSE
    )
  })
  df_w <- rbindlist(res_list)
  df_w$wilcox_fdr <- p.adjust(df_w$wilcox_p, method = "BH")

  # 2) LinDA compositional differential abundance (feature.dat: taxa x samples)
  feat_mat <- t(mat_filt)
  linda_fit <- tryCatch({
    linda(
      feature.dat = feat_mat,
      meta.dat = meta,
      formula = "~ Group",
      feature.dat.type = "proportion",
      prev.filter = 0.10,
      is.winsor = TRUE,
      zero.handling = "pseudo-count",
      p.adj.method = "BH"
    )
  }, error = function(e) {
    cat("  LinDA warning:", e$message, "\n")
    NULL
  })

  if (!is.null(linda_fit)) {
    l_out <- linda_fit$output[[1]]
    l_out$genus <- rownames(l_out)
    df_w <- merge(df_w, l_out[, c("genus", "log2FoldChange", "pvalue", "padj")], by = "genus", all.x = TRUE)
    colnames(df_w)[colnames(df_w) %in% c("log2FoldChange", "pvalue", "padj")] <- c("linda_log2FC", "linda_p", "linda_fdr")
  } else {
    df_w$linda_log2FC <- NA
    df_w$linda_p <- NA
    df_w$linda_fdr <- NA
  }

  all_diff[[proj]] <- df_w
  sig_w <- df_w[df_w$wilcox_fdr < 0.05, ]
  sig_l <- df_w[!is.na(df_w$linda_fdr) & df_w$linda_fdr < 0.05, ]
  cat(sprintf("  Significant genera: Wilcoxon FDR<0.05 = %d | LinDA FDR<0.05 = %d\n", nrow(sig_w), nrow(sig_l)))
}

df_alpha <- rbindlist(alpha_summary)
fwrite(df_alpha, file.path(res_dir, "GMrepo_alpha_beta_summary.csv"))
cat("\n=== Alpha & Beta Diversity Summary ===\n")
print(df_alpha)

df_all <- rbindlist(all_diff, fill = TRUE)
fwrite(df_all, file.path(res_dir, "GMrepo_genus_differential_by_cohort.csv"))

# Cross-cohort consistency summary
df_sig <- df_all[df_all$wilcox_p < 0.05 | (!is.na(df_all$linda_p) & df_all$linda_p < 0.05), ]
cons <- df_all[, .(
  n_cohorts_present = .N,
  n_wilcox_p05 = sum(wilcox_p < 0.05, na.rm = TRUE),
  n_wilcox_fdr05 = sum(wilcox_fdr < 0.05, na.rm = TRUE),
  n_linda_p05 = sum(linda_p < 0.05, na.rm = TRUE),
  n_linda_fdr05 = sum(linda_fdr < 0.05, na.rm = TRUE),
  mean_log2FC = mean(log2FC, na.rm = TRUE),
  dir_consistent = (all(log2FC > 0) || all(log2FC < 0)),
  cohorts_up = paste(project[log2FC > 0 & (wilcox_p < 0.05 | (!is.na(linda_p) & linda_p < 0.05))], collapse = ";"),
  cohorts_down = paste(project[log2FC < 0 & (wilcox_p < 0.05 | (!is.na(linda_p) & linda_p < 0.05))], collapse = ";")
), by = genus]
cons <- cons[order(-n_wilcox_fdr05, -n_linda_fdr05, -n_wilcox_p05)]
fwrite(cons, file.path(res_dir, "GMrepo_genus_cross_cohort_consistency.csv"))
cat("\n=== Top 25 Cross-Cohort Consistent Differential Genera ===\n")
print(head(cons, 25))
