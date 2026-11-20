project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(survival)
  library(ggrepel)
  library(ragg)
  library(systemfonts)
  library(grid)
})

arial_dir <- Sys.getenv("ARIAL_FONT_DIR")
if (nzchar(arial_dir)) {
  register_font(
    name = "Arial",
    plain = file.path(arial_dir, "arial.ttf"),
    bold = file.path(arial_dir, "arialbd.ttf"),
    italic = file.path(arial_dir, "ariali.ttf"),
    bolditalic = file.path(arial_dir, "arialbi.ttf")
  )
}

res_dir   <- file.path(project_root, "02_results")
fig_dir   <- file.path(project_root, "03_figures")
cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

PAL <- list(
  low   = "#2B7A9E",
  high  = "#C84B31",
  emer  = "#1F8A70",
  navy  = "#1E3A5F",
  slate = "#64748B",
  ink   = "#1A202C",
  muted = "#4A5568",
  line  = "#2D3748",
  div_lo = "#1D5B9B",
  div_mi = "#F8FAFC",
  div_hi = "#C0392B"
)

theme_top <- function(base_size = 8.2) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", color = PAL$ink),
      axis.line = element_line(color = PAL$line, linewidth = 0.45),
      axis.ticks = element_line(color = PAL$line, linewidth = 0.4),
      axis.ticks.length = unit(2.4, "pt"),
      axis.title = element_text(face = "bold", size = base_size, color = PAL$ink),
      axis.text = element_text(size = base_size - 0.8, color = "#2D3748"),
      plot.title = element_text(face = "bold", size = base_size + 1.2, hjust = 0, margin = margin(b = 4)),
      legend.background = element_blank(),
      legend.key = element_blank(),
      legend.title = element_text(face = "bold", size = base_size - 0.8),
      legend.text = element_text(size = base_size - 1.0),
      strip.background = element_blank(),
      strip.text = element_text(face = "bold", size = base_size - 0.2, color = PAL$ink),
      plot.margin = margin(5, 5, 4, 5)
    )
}

# ==============================================================================
# PANEL A: Compact Model Benchmark (UpSet Gene Dot-Matrix + C/AUC Heatmap)
# ==============================================================================
bench <- fread(file.path(res_dir, "Signature_Candidate_Selection_All5Cohorts.csv"))
top_models <- c(
  "C16_IfngIdo1Cxcl10Eno1_Ridge",
  "C16_IfngIdo1Cxcl10Eno1_UniBeta",
  "C1_IfngIdo1Cxcl10Tgfb1_UniBeta",
  "C18_IfngIdo1Cxcl10Eno1Tgfb1_UniBeta",
  "C22_IfngIdo1Cxcl10Tgfb1Tcf7l2_UniBeta",
  "C25_IfngIdo1Cxcl10Eno1Pemt_UniBeta",
  "C10_IfngIdo1Cxcl10VegfaTgfb1_UniBeta",
  "C2_IfngIdo1Cxcl10Tgfbr1_UniBeta",
  "C4_IfngIdo1Cxcl10Tgfb1Cldn7_UniBeta",
  "C8_IfngIdo1Cxcl10Tgfb1Magi3_UniBeta"
)
bsub <- bench[model %in% top_models]
model_names <- c(
  "C16_IfngIdo1Cxcl10Eno1_Ridge"          = "GMSig-4 (Ridge)*",
  "C16_IfngIdo1Cxcl10Eno1_UniBeta"        = "GMSig-4 (UniBeta)",
  "C1_IfngIdo1Cxcl10Tgfb1_UniBeta"        = "4-Gene (TGFB1)",
  "C18_IfngIdo1Cxcl10Eno1Tgfb1_UniBeta"   = "5-Gene (ENO1+TGFB1)",
  "C22_IfngIdo1Cxcl10Tgfb1Tcf7l2_UniBeta" = "5-Gene (+TCF7L2)",
  "C25_IfngIdo1Cxcl10Eno1Pemt_UniBeta"    = "5-Gene (+PEMT)",
  "C10_IfngIdo1Cxcl10VegfaTgfb1_UniBeta" = "5-Gene (+VEGFA)",
  "C2_IfngIdo1Cxcl10Tgfbr1_UniBeta"       = "4-Gene (TGFBR1)",
  "C4_IfngIdo1Cxcl10Tgfb1Cldn7_UniBeta"  = "5-Gene (+CLDN7)",
  "C8_IfngIdo1Cxcl10Tgfb1Magi3_UniBeta"  = "5-Gene (+MAGI3)"
)
mod_ord <- rev(unname(model_names[top_models]))
bsub[, label := factor(model_names[model], levels = mod_ord)]

all_genes <- c("IFNG", "IDO1", "CXCL10", "ENO1", "TGFB1", "TCF7L2", "PEMT", "VEGFA", "TGFBR1", "CLDN7", "MAGI3")
gene_membership <- list(
  "GMSig-4 (Ridge)*"     = c("IFNG", "IDO1", "CXCL10", "ENO1"),
  "GMSig-4 (UniBeta)"    = c("IFNG", "IDO1", "CXCL10", "ENO1"),
  "4-Gene (TGFB1)"       = c("IFNG", "IDO1", "CXCL10", "TGFB1"),
  "5-Gene (ENO1+TGFB1)"  = c("IFNG", "IDO1", "CXCL10", "ENO1", "TGFB1"),
  "5-Gene (+TCF7L2)"     = c("IFNG", "IDO1", "CXCL10", "TGFB1", "TCF7L2"),
  "5-Gene (+PEMT)"       = c("IFNG", "IDO1", "CXCL10", "ENO1", "PEMT"),
  "5-Gene (+VEGFA)"      = c("IFNG", "IDO1", "CXCL10", "VEGFA", "TGFB1"),
  "4-Gene (TGFBR1)"      = c("IFNG", "IDO1", "CXCL10", "TGFBR1"),
  "5-Gene (+CLDN7)"      = c("IFNG", "IDO1", "CXCL10", "TGFB1", "CLDN7"),
  "5-Gene (+MAGI3)"      = c("IFNG", "IDO1", "CXCL10", "TGFB1", "MAGI3")
)
up_dt <- CJ(label = factor(mod_ord, levels = mod_ord), gene = factor(all_genes, levels = all_genes))
up_dt[, present := mapply(function(l, g) g %in% gene_membership[[as.character(l)]], label, gene)]
up_dt[, gx := as.numeric(gene)]
up_lines <- up_dt[present == TRUE, .(xmin = min(gx), xmax = max(gx)), by = label]

p3a_upset <- ggplot(up_dt, aes(x = gene, y = label)) +
  geom_tile(fill = "#F8FAFC", color = "white", linewidth = 0.5) +
  geom_segment(data = up_lines, aes(x = xmin, xend = xmax, y = label, yend = label),
               inherit.aes = FALSE, color = PAL$navy, linewidth = 0.75) +
  geom_point(data = up_dt[present == FALSE], size = 1.9, color = "#CBD5E0") +
  geom_point(data = up_dt[present == TRUE], size = 2.6,
             aes(color = ifelse(grepl("GMSig-4", label), PAL$high, PAL$navy))) +
  scale_color_identity() +
  scale_x_discrete(labels = parse(text = sprintf("italic('%s')", all_genes))) +
  labs(x = NULL, y = NULL) +
  theme_top(8.0) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(angle = 38, hjust = 1, size = 6.9, face = "bold"),
    axis.text.y = element_text(size = 7.3, face = "bold"),
    plot.margin = margin(6, 2, 4, 4)
  )

mlong <- melt(
  bsub[, .(label,
           `TCGA\nOS (C)`       = tcga_C,
           `GSE58812\nOS (C)`   = g58_OS_C,
           `GSE96058\nExt (C)`  = g960e_C,
           `Pembro\nTNBC`       = pembro_TN_AUC,
           `Pembro\nAll`        = pembro_All_AUC,
           `Durva\nAll`         = durva_All_AUC,
           `Durva\nTNBC`        = durva_TN_AUC,
           `Paclitaxel\nCtrl`   = paclitaxel_TN_AUC)],
  id.vars = "label", variable.name = "Metric", value.name = "Score"
)
mlong[, txt_col := ifelse(!is.na(Score) & (Score >= 0.76 | Score <= 0.50), "white", PAL$ink)]

p3a_heat <- ggplot(mlong, aes(x = Metric, y = label, fill = Score)) +
  geom_tile(color = "white", linewidth = 0.65) +
  geom_tile(data = mlong[label == "GMSig-4 (Ridge)*"], fill = NA, color = PAL$ink, linewidth = 0.7) +
  geom_text(aes(label = ifelse(is.na(Score), "NA", sprintf("%.2f", Score)), color = txt_col),
            size = 2.4, fontface = "bold", family = "Arial") +
  scale_color_identity() +
  scale_fill_gradient2(low = PAL$div_lo, mid = PAL$div_mi, high = PAL$div_hi, midpoint = 0.64,
                       limits = c(0.45, 0.85), na.value = "#EDF2F7", name = "C / AUC") +
  labs(x = NULL, y = NULL) +
  theme_top(8.0) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_blank(), axis.text.x = element_text(size = 6.9, face = "bold"),
    legend.position = "right", plot.margin = margin(6, 4, 4, 1)
  ) +
  guides(fill = guide_colorbar(barwidth = unit(6, "pt"), barheight = unit(72, "pt")))

p3a <- (p3a_upset | p3a_heat) + plot_layout(widths = c(1.05, 1.45))

# ==============================================================================
# PANEL B: Locked GMSig-4 Ridge Coefficients
# ==============================================================================
m_obj <- readRDS(file.path(res_dir, "GMSig4_Model_Objects.rds"))
w_vec <- m_obj$weights
w_dt <- data.table(gene = names(w_vec), beta = as.numeric(w_vec))
w_dt[, Dir := ifelse(beta > 0, "Risk (+)", "Protective (\u2212)")]
w_dt <- w_dt[order(beta)]
w_dt[, gene_f := factor(gene, levels = gene)]

p3b <- ggplot(w_dt, aes(x = beta, y = gene_f, color = Dir, fill = Dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = beta, yend = gene_f), linewidth = 1.1) +
  geom_point(shape = 21, size = 3.8, color = "white", stroke = 0.6) +
  geom_text(aes(label = sprintf("%+.4f", beta), hjust = ifelse(beta >= 0, -0.25, 1.25)),
            size = 2.55, fontface = "bold", family = "Arial", color = PAL$ink) +
  annotate("label", x = 0.14, y = 1.25, label = "TCGA-TNBC 5-fold CV\nEPV = 29 / 4 = 7.25",
           size = 2.3, fontface = "bold", family = "Arial", color = PAL$navy,
           fill = "#F8FAFC", linewidth = 0.3) +
  scale_color_manual(values = c("Protective (\u2212)" = PAL$low, "Risk (+)" = PAL$high), guide = "none") +
  scale_fill_manual(values = c("Protective (\u2212)" = PAL$low, "Risk (+)" = PAL$high), name = NULL) +
  scale_x_continuous(limits = c(-0.37, 0.29), breaks = c(-0.3, -0.15, 0, 0.15)) +
  scale_y_discrete(labels = parse(text = sprintf("italic('%s')", levels(w_dt$gene_f)))) +
  labs(x = "Ridge Cox Weight (\u03B2)", y = NULL) +
  theme_top(8.1) +
  theme(axis.text.y = element_text(size = 8.4, face = "bold"),
        legend.position = "bottom", legend.margin = margin(t = -3))

# ==============================================================================
# PANELS C-F: Kaplan-Meier Curves + Aligned Number-at-Risk Tables (Zero Dup Titles)
# ==============================================================================
tcga_df <- fread(file.path(res_dir, "GMSig4_PatientScores_TCGA_TNBC.csv"), data.table = FALSE)
g58_df  <- fread(file.path(res_dir, "GMSig4_PatientScores_GSE58812.csv"), data.table = FALSE)
g960_df <- fread(file.path(res_dir, "GMSig4_PatientScores_GSE96058_TNBC_Ext.csv"), data.table = FALSE)
ssum    <- fread(file.path(res_dir, "GMSig4_Prognostic_Validation_Summary.csv"))

fmt_km_lbl <- function(row_dt, extra = NULL) {
  p_str <- ifelse(row_dt$KM_logrank_p < 1e-3,
                  formatC(row_dt$KM_logrank_p, format = "e", digits = 2),
                  sprintf("%.4f", row_dt$KM_logrank_p))
  base <- sprintf("n=%d (%d events) | C=%.3f\nHR=%.2f (%.2f\u2013%.2f)\nLog-rank P=%s",
                  row_dt$n, row_dt$events, row_dt$C_index,
                  row_dt$HighLow_HR, row_dt$HighLow_LCI, row_dt$HighLow_UCI, p_str)
  if (!is.null(extra)) paste0(base, "\n", extra) else base
}

make_km_with_risk <- function(df, tcol, ecol, stats_lbl, show_legend = FALSE) {
  d <- df[!is.na(df[[tcol]]) & !is.na(df[[ecol]]) & df[[tcol]] > 0, ]
  d$Grp <- factor(ifelse(d$GMSig4 > median(d$GMSig4), "High Risk", "Low Risk"), levels = c("Low Risk", "High Risk"))
  sf <- survfit(as.formula(paste0("Surv(", tcol, ", ", ecol, ") ~ Grp")), data = d)
  km_df <- data.frame(
    time = sf$time, surv = sf$surv, lower = sf$lower, upper = sf$upper,
    strata = factor(rep(c("Low Risk", "High Risk"), sf$strata), levels = c("Low Risk", "High Risk"))
  )
  km_df <- rbind(
    data.frame(time = c(0, 0), surv = c(1, 1), lower = c(1, 1), upper = c(1, 1),
               strata = factor(c("Low Risk", "High Risk"), levels = c("Low Risk", "High Risk"))),
    km_df
  )
  max_t <- ceiling(max(d[[tcol]], na.rm = TRUE))
  t_breaks <- pretty(c(0, max_t), n = 4)
  t_breaks <- t_breaks[t_breaks <= max_t + 1]
  sm <- summary(sf, times = t_breaks, extend = TRUE)
  rt_df <- data.frame(
    time = sm$time,
    strata = factor(gsub("Grp=", "", as.character(sm$strata)), levels = c("High Risk", "Low Risk")),
    n.risk = sm$n.risk
  )

  p_curve <- ggplot(km_df, aes(x = time, y = surv, color = strata, fill = strata)) +
    geom_ribbon(aes(ymin = lower, ymax = upper), alpha = 0.13, color = NA) +
    geom_step(linewidth = 0.85) +
    scale_color_manual(values = c("Low Risk" = PAL$low, "High Risk" = PAL$high), name = NULL) +
    scale_fill_manual(values = c("Low Risk" = PAL$low, "High Risk" = PAL$high), name = NULL) +
    scale_x_continuous(breaks = t_breaks, limits = c(0, max(t_breaks))) +
    scale_y_continuous(limits = c(0, 1.03), breaks = seq(0, 1, 0.25)) +
    annotate("text", x = max(t_breaks) * 0.03, y = 0.04, hjust = 0, vjust = 0,
             label = stats_lbl, size = 2.25, fontface = "bold", family = "Arial", color = PAL$ink) +
    labs(x = NULL, y = "Survival Probability") +
    theme_top(8.0) +
    theme(
      legend.position = if (show_legend) c(0.75, 0.88) else "none",
      plot.margin = margin(10, 5, 2, 5)
    )

  p_risk <- ggplot(rt_df, aes(x = time, y = strata, label = n.risk, color = strata)) +
    geom_text(size = 2.25, fontface = "bold", family = "Arial") +
    scale_color_manual(values = c("Low Risk" = PAL$low, "High Risk" = PAL$high), guide = "none") +
    scale_x_continuous(breaks = t_breaks, limits = c(0, max(t_breaks))) +
    labs(x = "Time (Years)", y = NULL) +
    theme_top(7.6) +
    theme(
      axis.line.y = element_blank(), axis.ticks.y = element_blank(),
      axis.text.y = element_text(size = 6.7, face = "bold", color = c(PAL$high, PAL$low)),
      plot.margin = margin(0, 5, 4, 5)
    )

  (p_curve / p_risk) + plot_layout(heights = c(1.0, 0.26))
}

p3c <- make_km_with_risk(tcga_df, "OS_yr", "OS", fmt_km_lbl(ssum[cohort == "TCGA_TNBC_OS"]), show_legend = TRUE)
p3d <- make_km_with_risk(g58_df,  "OS_yr", "OS", fmt_km_lbl(ssum[cohort == "GSE58812_TNBC_OS"]))
p3e <- make_km_with_risk(g58_df,  "MFS_yr","MFS",fmt_km_lbl(ssum[cohort == "GSE58812_TNBC_MFS"]))
p3f <- make_km_with_risk(g960_df, "OS_yr", "OS", fmt_km_lbl(ssum[cohort == "GSE96058_TNBC_Ext_OS"], "Strict (n=143): P=0.0051"))

# ==============================================================================
# PANEL G: Time-Dependent ROC Trajectory w/ ggrepel Zero Overlap
# ==============================================================================
roc_dt <- melt(
  ssum[cohort %in% c("TCGA_TNBC_OS", "GSE58812_TNBC_OS", "GSE96058_TNBC_Strict_OS", "GSE96058_TNBC_Ext_OS"),
       .(cohort, `1-Year` = AUC_1yr, `3-Year` = AUC_3yr, `5-Year` = AUC_5yr)],
  id.vars = "cohort", variable.name = "Horizon", value.name = "AUC"
)
roc_dt[, Cohort := factor(c(
  "TCGA_TNBC_OS"            = "TCGA-TNBC (n=172)",
  "GSE58812_TNBC_OS"        = "GSE58812 (n=107)",
  "GSE96058_TNBC_Strict_OS" = "GSE96058 Strict (n=143)",
  "GSE96058_TNBC_Ext_OS"    = "GSE96058 Ext (n=242)"
)[cohort], levels = c("TCGA-TNBC (n=172)", "GSE58812 (n=107)", "GSE96058 Strict (n=143)", "GSE96058 Ext (n=242)"))]

p3g <- ggplot(roc_dt, aes(x = Horizon, y = AUC, group = Cohort, color = Cohort, shape = Cohort)) +
  geom_hline(yintercept = 0.5, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_line(linewidth = 0.9, alpha = 0.85) +
  geom_point(size = 2.9) +
  ggrepel::geom_text_repel(
    aes(label = sprintf("%.2f", AUC)),
    size = 2.25, fontface = "bold", family = "Arial",
    box.padding = 0.25, point.padding = 0.15, min.segment.length = 0.2,
    segment.size = 0.3, show.legend = FALSE, seed = 42
  ) +
  scale_color_manual(values = c(PAL$high, PAL$emer, PAL$low, PAL$navy), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15, 18), name = NULL) +
  scale_y_continuous(limits = c(0.45, 0.86), breaks = seq(0.5, 0.8, 0.1)) +
  labs(x = "Prediction Horizon", y = "Time-Dependent AUC") +
  theme_top(8.1) +
  theme(legend.position = "bottom", legend.margin = margin(t = -3),
        plot.margin = margin(10, 5, 4, 5)) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

# ==============================================================================
# PANEL H: Multivariable Independent Prognosis Forest Table (Split Forest + Table)
# ==============================================================================
mv_dt <- fread(file.path(res_dir, "GMSig4_Multivariable_Cox_Results.csv"))
mv_dt[, cohort_clean := factor(
  ifelse(cohort == "TCGA_TNBC", "TCGA-TNBC (n=168, C=0.866)", "GSE96058-TNBC (n=234, C=0.785)"),
  levels = c("TCGA-TNBC (n=168, C=0.866)", "GSE96058-TNBC (n=234, C=0.785)")
)]
mv_dt[, var_only := c("GMSig4" = "GMSig-4 Score", "Age" = "Age (per yr)",
                      "Stage_III_IV" = "Stage III/IV", "Node_Pos" = "Lymph Node+",
                      "Tumor_mm" = "Tumor Size (mm)")[variable]]
mv_dt[, var_f := reorder(paste(cohort, var_only, sep = "__"), HR)]
mv_dt[, tbl_or := sprintf("%.2f (%.2f\u2013%.2f)", HR, LCI, UCI)]
mv_dt[, tbl_p  := ifelse(pval < 0.001, formatC(pval, format = "e", digits = 1), sprintf("%.3f", pval))]
mv_dt[, is_sig := variable == "GMSig4"]

p3h_forest <- ggplot(mv_dt, aes(x = HR, y = var_f, color = is_sig)) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_errorbar(aes(xmin = LCI, xmax = UCI), orientation = "y", width = 0.22, linewidth = 0.7) +
  geom_point(size = 2.6) +
  facet_wrap(~cohort_clean, scales = "free_y", ncol = 1) +
  scale_x_log10(breaks = c(0.5, 1, 2, 5, 15), limits = c(0.42, 20)) +
  scale_y_discrete(labels = function(z) sub("^.*__", "", z)) +
  scale_color_manual(values = c("FALSE" = PAL$navy, "TRUE" = PAL$high), guide = "none") +
  labs(x = "Hazard Ratio (log scale)", y = NULL) +
  theme_top(8.1) +
  theme(
    strip.background = element_rect(fill = "#EDF2F7", color = NA),
    strip.text = element_text(face = "bold", size = 7.6, color = PAL$navy),
    plot.margin = margin(14, 1, 4, 5)
  )

p3h_tbl <- ggplot(mv_dt, aes(y = var_f)) +
  geom_text(aes(x = 0.0, label = tbl_or, fontface = ifelse(is_sig, "bold", "plain")),
            hjust = 0, size = 2.25, family = "Arial", color = PAL$ink) +
  geom_text(aes(x = 1.05, label = tbl_p, fontface = ifelse(is_sig, "bold", "plain"),
                color = ifelse(is_sig, PAL$high, PAL$ink)),
            hjust = 1, size = 2.25, family = "Arial") +
  facet_wrap(~cohort_clean, scales = "free_y", ncol = 1) +
  scale_color_identity() +
  scale_x_continuous(limits = c(-0.02, 1.10), breaks = c(0.38, 0.96), labels = c("HR (95% CI)", "P")) +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 7.1, color = PAL$muted),
    strip.background = element_rect(fill = "#EDF2F7", color = NA),
    strip.text = element_text(face = "bold", size = 7.6, color = "#EDF2F7"),
    plot.margin = margin(14, 4, 4, 0)
  )

p3h <- (p3h_forest | p3h_tbl) + plot_layout(widths = c(1.22, 0.95))

# ==============================================================================
# PANEL I: Bootstrap 1/3/5-Year Calibration Curves w/ 95% CI Error Bars
# ==============================================================================
cal_dt <- fread(file.path(res_dir, "GMSig4_Nomogram_Calibration_1_3_5yr.csv"))
cal_dt[, Horizon := factor(paste0(time_yr, "-Year OS"), levels = c("1-Year OS", "3-Year OS", "5-Year OS"))]
cal_dt[, Cohort := factor(
  ifelse(cohort == "TCGA_TNBC", "TCGA-TNBC (C=0.866, LRT P=0.0081)", "GSE96058 (C=0.785, LRT P=0.0128)"),
  levels = c("TCGA-TNBC (C=0.866, LRT P=0.0081)", "GSE96058 (C=0.785, LRT P=0.0128)")
)]
cal_dt[, lci := pmax(0, KM - 1.96 * std.err)]
cal_dt[, uci := pmin(1, KM + 1.96 * std.err)]

p3i <- ggplot(cal_dt, aes(x = mean.predicted, y = KM, color = Horizon)) +
  geom_abline(slope = 1, intercept = 0, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_errorbar(aes(ymin = lci, ymax = uci), width = 0.015, linewidth = 0.45, alpha = 0.55) +
  geom_line(linewidth = 0.8) +
  geom_point(size = 2.0) +
  facet_wrap(~Cohort, ncol = 1) +
  scale_color_manual(values = c("1-Year OS" = PAL$low, "3-Year OS" = PAL$emer, "5-Year OS" = PAL$high), name = NULL) +
  coord_cartesian(xlim = c(0.38, 1.02), ylim = c(0.25, 1.04)) +
  labs(x = "Predicted Survival Probability", y = "Observed KM Survival") +
  theme_top(8.1) +
  theme(
    strip.background = element_rect(fill = "#EDF2F7", color = NA),
    strip.text = element_text(face = "bold", size = 7.6, color = PAL$navy),
    legend.position = "bottom", legend.margin = margin(t = -3),
    plot.margin = margin(14, 5, 4, 5)
  )

# ==============================================================================
# Clean Composite Assembly (Panels A-I)
# ==============================================================================
wrap_panel <- function(plot_obj, tag_letter, title_txt) {
  wrap_elements(full = plot_obj) +
    labs(tag = tag_letter, title = title_txt) +
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 14, color = "#111111"),
      plot.tag.position = c(0.01, 0.986),
      plot.title = element_text(family = "Arial", face = "bold", size = 9.5, color = PAL$ink,
                                hjust = 0, margin = margin(l = 22, t = 2, b = 7)),
      plot.margin = margin(4, 4, 4, 4)
    )
}

w3a <- wrap_panel(p3a, "A", "Multi-Cohort Candidate Signature Benchmark (3\u20135 Genes)")
w3b <- wrap_panel(p3b, "B", "Locked GMSig-4 Ridge Weights")
w3c <- wrap_panel(p3c, "C", "TCGA-TNBC OS")
w3d <- wrap_panel(p3d, "D", "GSE58812 OS")
w3e <- wrap_panel(p3e, "E", "GSE58812 MFS")
w3f <- wrap_panel(p3f, "F", "GSE96058 Extended OS")
w3g <- wrap_panel(p3g, "G", "Time-Dependent ROC AUCs")
w3h <- wrap_panel(p3h, "H", "Multivariable Independent Prognosis")
w3i <- wrap_panel(p3i, "I", "Nomogram 1/3/5-Yr Calibration")

row1 <- (w3a | w3b) + plot_layout(widths = c(2.05, 0.85))
row2 <- (w3c | w3d | w3e | w3f)
row3 <- (w3g | w3h | w3i) + plot_layout(widths = c(0.92, 1.32, 0.92))

fig3 <- row1 / row2 / row3 + plot_layout(heights = c(1.02, 1.15, 1.18))

agg_png(file.path(fig_dir, "Figure3.png"), width = 16.4, height = 15.6, units = "in", res = 300, background = "white")
print(fig3)
dev.off()

cairo_pdf(file.path(fig_dir, "Figure3.pdf"), width = 16.4, height = 15.6)
print(fig3)
dev.off()

cat("SAVED Figure3.png and Figure3.pdf successfully!\n")
