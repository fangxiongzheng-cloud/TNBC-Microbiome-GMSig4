#!/usr/bin/env Rscript
# Publication-Grade Figure 6: Strictly Independent External Validation,
# Paired Bootstrap Comparator Benchmarks & Optimism-Corrected Clinical Incremental Value.

.libPaths(c("/root/miniconda3/envs/bioinfo-r-core/lib/R/library", "/root/R/tnbc_microbiome_lib", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(cowplot)
})

root    <- "/mnt/d/Bioinformatics/TNBC_Microbiome"
val_dir <- file.path(root, "02_results/JTM_upgrade/independent_validation")
fig_dir <- file.path(root, "03_figures")

PAL <- list(
  high  = "#C0392B",
  low   = "#1B6B93",
  navy  = "#1B365D",
  emer  = "#2E8B57",
  amber = "#D9822B",
  muted = "#64748B",
  ink   = "#1E293B",
  grid  = "#E2E8F0"
)

theme_pub <- function(base_size = 9.0) {
  theme_classic(base_size = base_size) +
    theme(
      text = element_text(color = PAL$ink),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(size = base_size - 0.8, color = PAL$ink),
      axis.line = element_line(color = PAL$ink, linewidth = 0.45),
      axis.ticks = element_line(color = PAL$ink, linewidth = 0.4),
      plot.title = element_text(face = "bold", size = base_size + 0.7, color = PAL$navy, hjust = 0),
      plot.subtitle = element_text(size = base_size - 1.1, color = PAL$muted, hjust = 0),
      plot.title.position = "plot",
      legend.background = element_blank(),
      legend.key = element_blank(),
      strip.background = element_rect(fill = "#EDF2F7", color = NA),
      strip.text = element_text(face = "bold", size = base_size - 0.8, color = PAL$navy),
      plot.margin = margin(6, 8, 6, 8)
    )
}

make_km_panel <- function(df, time_col, event_col, title_txt, sub_txt, annot_txt, xmax = NULL) {
  d <- copy(df[!is.na(get(time_col)) & !is.na(get(event_col)) & get(time_col) > 0])
  d[, grp := factor(RiskGroup, levels = c("Low", "High"))]
  sf <- survfit(Surv(d[[time_col]], d[[event_col]]) ~ grp, data = d)
  s_dt <- rbindlist(lapply(seq_along(sf$strata), function(i) {
    idx <- if (i == 1) seq_len(sf$strata[1]) else (sum(sf$strata[seq_len(i - 1)]) + 1):sum(sf$strata[seq_len(i)])
    lbl <- sub("^grp=", "", names(sf$strata)[i])
    data.table(
      time = c(0, sf$time[idx]),
      surv = c(1, sf$surv[idx]),
      lower = c(1, sf$lower[idx]),
      upper = c(1, sf$upper[idx]),
      grp = lbl
    )
  }))
  s_dt[, grp := factor(grp, levels = c("Low", "High"))]
  if (!is.null(xmax)) s_dt <- s_dt[time <= xmax]
  n_low  <- sum(d$grp == "Low")
  n_high <- sum(d$grp == "High")
  lbl_map <- c(
    "Low"  = sprintf("Low GMSig-4 (n=%d)", n_low),
    "High" = sprintf("High GMSig-4 (n=%d)", n_high)
  )
  ggplot(s_dt, aes(x = time, y = surv, color = grp, fill = grp)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.14, color = NA, stat = "identity") +
    geom_step(linewidth = 0.95) +
    scale_color_manual(values = c("Low" = PAL$low, "High" = PAL$high), labels = lbl_map, name = NULL) +
    scale_fill_manual(values = c("Low" = PAL$low, "High" = PAL$high), guide = "none") +
    scale_y_continuous(limits = c(0, 1.02), breaks = seq(0, 1, 0.25), expand = c(0.01, 0)) +
    scale_x_continuous(expand = c(0.02, 0)) +
    annotate("label", x = if (!is.null(xmax)) xmax * 0.03 else max(s_dt$time) * 0.03,
             y = 0.16, label = annot_txt, hjust = 0, vjust = 0.5, size = 2.65,
             fill = alpha("white", 0.92), color = PAL$ink, linewidth = 0.25, lineheight = 1.15) +
    labs(title = title_txt, subtitle = sub_txt, x = "Follow-Up Time (Years)", y = "Survival Probability") +
    theme_pub(8.8) +
    theme(legend.position = c(0.68, 0.90), legend.text = element_text(size = 7.4, face = "bold"))
}

# Load Stage B outputs
sum_dt  <- fread(file.path(val_dir, "Independent_Validation_Summary.csv"))
comp_dt <- fread(file.path(val_dir, "Independent_Comparator_Benchmarks.csv"))
inc_dt  <- fread(file.path(val_dir, "Independent_Clinical_Incremental_Value.csv"))
mb_df   <- fread(file.path(val_dir, "PatientScores_METABRIC_TNBC_OS.csv"))
g21_df  <- fread(file.path(val_dir, "PatientScores_GSE21653_TNBC_DFS.csv"))
g25_df  <- fread(file.path(val_dir, "PatientScores_GSE25066_TNBC_DRFS.csv"))

# Panels A, B, C: KM Curves
r_mb <- sum_dt[cohort_endpoint == "METABRIC_TNBC_OS"]
r_mb_dss <- sum_dt[cohort_endpoint == "METABRIC_TNBC_DSS"]
r_mb_rfs <- sum_dt[cohort_endpoint == "METABRIC_TNBC_RFS"]
pA <- make_km_panel(
  mb_df, "OS_yr", "OS_event",
  "A   METABRIC Primary TNBC OS (n = 319, 168 events)",
  "Strictly independent cohort (0% prior exposure; rev-KM 15.2 yr)",
  sprintf("Per-SD HR = %.2f (%.2f-%.2f), Cox P = %.4f\nHigh/Low HR = %.2f, Log-rank P = %.4f\nC = %.3f | DSS P = %.4f | RFS P = %.4f",
          r_mb$HR_per_cohort_SD, r_mb$HR_SD_LCI, r_mb$HR_SD_UCI, r_mb$Cox_p,
          r_mb$HighLow_HR, r_mb$KM_logrank_p, r_mb$C_index, r_mb_dss$Cox_p, r_mb_rfs$Cox_p),
  xmax = 18
)

r_21 <- sum_dt[cohort_endpoint == "GSE21653_TNBC_DFS"]
pB <- make_km_panel(
  g21_df, "DFS_yr", "DFS_event",
  "B   GSE21653 Strict IHC TNBC DFS (n = 83, 25 events)",
  "Independent IPC Marseille cohort (0% overlap with GSE58812)",
  sprintf("Per-SD HR = %.2f (%.2f-%.2f), Cox P = %.1e\nHigh/Low HR = %.2f, Log-rank P = %.4f\nC = %.3f (%.3f-%.3f) | 3-yr AUC = %.3f",
          r_21$HR_per_cohort_SD, r_21$HR_SD_LCI, r_21$HR_SD_UCI, r_21$Cox_p,
          r_21$HighLow_HR, r_21$KM_logrank_p, r_21$C_index, r_21$C_LCI, r_21$C_UCI, r_21$AUC_3yr),
  xmax = 10
)

r_25 <- sum_dt[cohort_endpoint == "GSE25066_TNBC_DRFS"]
pC <- make_km_panel(
  g25_df, "DRFS_yr", "DRFS_event",
  "C   GSE25066 Neoadjuvant Chemo DRFS (n = 178, 64 events)",
  "Honest null in taxane-anthracycline cohort (chemo pCR AUC = 0.442)",
  sprintf("Per-SD HR = %.2f (%.2f-%.2f), Cox P = %.3f\nHigh/Low HR = %.2f, Log-rank P = %.3f\nC = %.3f (all immune comparators P > 0.16)",
          r_25$HR_per_cohort_SD, r_25$HR_SD_LCI, r_25$HR_SD_UCI, r_25$Cox_p,
          r_25$HighLow_HR, r_25$KM_logrank_p, r_25$C_index),
  xmax = 6.5
)

# Panel D: Paired Bootstrap Delta C-Index Forest Plot
comp_sub <- copy(comp_dt[cohort_endpoint %in% c("METABRIC_TNBC_OS", "METABRIC_TNBC_DSS", "METABRIC_TNBC_RFS", "GSE21653_TNBC_DFS")])
comp_sub[, ep_lbl := factor(
  c("METABRIC_TNBC_OS" = "METABRIC OS (n=319)",
    "METABRIC_TNBC_DSS" = "METABRIC DSS (n=319)",
    "METABRIC_TNBC_RFS" = "METABRIC RFS (n=319)",
    "GSE21653_TNBC_DFS" = "GSE21653 DFS (n=83)")[cohort_endpoint],
  levels = c("METABRIC OS (n=319)", "METABRIC DSS (n=319)", "METABRIC RFS (n=319)", "GSE21653 DFS (n=83)")
)]
comp_map <- c(
  "STAT1_single"     = "vs STAT1 (single gene)",
  "PDCD1_single"     = "vs PDCD1 (single gene)",
  "CD274_single"     = "vs CD274 (single gene)",
  "Immune3_frozen"   = "vs Immune-3 (no ENO1)",
  "ENO1_only"        = "vs ENO1 only",
  "GEP18_unweighted" = "vs 18-gene GEP",
  "CYT_GZMA_PRF1"    = "vs CYT (GZMA+PRF1)"
)
comp_sub[, comp_lbl := factor(comp_map[comparator], levels = rev(unname(comp_map)))]
comp_sub[, sig_tier := factor(
  ifelse(p_holm < 0.05, "Holm P < 0.05", ifelse(p_boot_two_sided < 0.05, "Nominal boot P < 0.05", "Parity (NS)")),
  levels = c("Holm P < 0.05", "Nominal boot P < 0.05", "Parity (NS)")
)]

pD <- ggplot(comp_sub, aes(x = Delta_C_GMSig4_minus_Comp, y = comp_lbl, color = sig_tier)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#64748B", linewidth = 0.45) +
  geom_errorbar(aes(xmin = Delta_C_LCI, xmax = Delta_C_UCI), orientation = "y", width = 0.25, linewidth = 0.65) +
  geom_point(size = 2.2) +
  facet_wrap(~ep_lbl, nrow = 1, scales = "free_x") +
  scale_color_manual(values = c("Holm P < 0.05" = PAL$high, "Nominal boot P < 0.05" = PAL$amber, "Parity (NS)" = PAL$navy), name = NULL) +
  labs(
    title = "D   Paired Bootstrap Discrimination (Delta C = GMSig-4 minus Comparator; B = 1,000)",
    subtitle = "Outperforms single-gene STAT1 & Immune-3 ablation (confirming ENO1 gain); parity with 18-gene GEP",
    x = "Paired Delta C-index (95% Bootstrap Percentile CI)", y = NULL
  ) +
  theme_pub(8.6) +
  theme(legend.position = "bottom", legend.margin = margin(t = -3),
        axis.text.y = element_text(face = "bold", size = 7.8),
        plot.margin = margin(6, 10, 6, 12))

# Panel E: Optimism-Corrected Clinical Incremental Value
inc_sub <- copy(inc_dt[cohort_model %in% c(
  "METABRIC_TNBC_OS_PrimaryAdj", "METABRIC_TNBC_DSS_PrimaryAdj",
  "METABRIC_TNBC_RFS_PrimaryAdj", "GSE25066_TNBC_DRFS_PrimaryAdj",
  "GSE21653_TNBC_DFS_Exploratory"
)])
inc_map <- c(
  "METABRIC_TNBC_OS_PrimaryAdj"   = "METABRIC OS\n(n=295, Age+Size+Node)",
  "METABRIC_TNBC_DSS_PrimaryAdj"  = "METABRIC DSS\n(n=295, Age+Size+Node)",
  "METABRIC_TNBC_RFS_PrimaryAdj"  = "METABRIC RFS\n(n=295, Age+Size+Node)",
  "GSE25066_TNBC_DRFS_PrimaryAdj" = "GSE25066 DRFS\n(n=176, Age+Node+Stage)",
  "GSE21653_TNBC_DFS_Exploratory" = "GSE21653 DFS*\n(n=82, EPV=8, Age+Node)"
)
inc_sub[, lbl := factor(inc_map[cohort_model], levels = rev(unname(inc_map)))]
inc_sub[, lrt_lbl := ifelse(LRT_p < 0.001, sprintf("LRT P=%.1e", LRT_p), sprintf("LRT P=%.4f", LRT_p))]

pE <- ggplot(inc_sub, aes(y = lbl)) +
  geom_segment(aes(x = C_clin_opt_corr, xend = C_full_opt_corr, yend = lbl), color = "#94A3B8", linewidth = 1.4) +
  geom_point(aes(x = C_clin_opt_corr, color = "Clinical Only"), size = 3.2) +
  geom_point(aes(x = C_full_opt_corr, color = "Clinical + GMSig-4"), size = 3.5) +
  geom_text(aes(x = pmax(C_clin_opt_corr, C_full_opt_corr) + 0.008,
                label = sprintf("%+.3f (%s)", Delta_C_opt_corr, lrt_lbl)),
            hjust = 0, size = 2.65, fontface = "bold", color = PAL$ink) +
  scale_color_manual(values = c("Clinical Only" = PAL$muted, "Clinical + GMSig-4" = PAL$high), name = NULL) +
  scale_x_continuous(limits = c(0.55, 0.93), breaks = seq(0.55, 0.90, 0.05)) +
  labs(
    title = "E   Clinical Incremental Value (500-Boot Optimism-Corrected)",
    subtitle = "Cohort-fitted clinical models + frozen GMSig-4 (*exploratory EPV < 10)",
    x = "Optimism-Corrected Harrell C-index", y = NULL
  ) +
  theme_pub(8.6) +
  theme(legend.position = "bottom", legend.margin = margin(t = -3),
        axis.text.y = element_text(size = 7.6, face = "bold"))

row1 <- plot_grid(pA, pB, pC, nrow = 1, align = "h")
row2 <- plot_grid(pD, pE, nrow = 1, rel_widths = c(1.42, 1.0))
fig6 <- plot_grid(row1, row2, ncol = 1, rel_heights = c(1.0, 1.05))

ggsave(file.path(fig_dir, "Figure6_Independent_Validation.png"), fig6, width = 15.2, height = 9.4, dpi = 300, bg = "white")
ggsave(file.path(fig_dir, "Figure6_Independent_Validation.pdf"), fig6, width = 15.2, height = 9.4, bg = "white")
cat("SAVED Figure6_Independent_Validation.png and .pdf successfully!\n")
