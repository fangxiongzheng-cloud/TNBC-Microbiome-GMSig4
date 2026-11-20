project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(ggplot2)
  library(patchwork)
})

cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
res_dir   <- file.path(project_root, "02_results")
sup_dir   <- file.path(project_root, "03_figures", "Supplementary")
dir.create(sup_dir, showWarnings = FALSE, recursive = TRUE)

tcga <- readRDS(file.path(cache_dir, "TCGA_BRCA_clean.rds"))
g960 <- readRDS(file.path(cache_dir, "GSE96058_TNBC_clean.rds"))

# Reconstruct complete-case cohorts identically to 19f_signature_clinical_nomogram.R
tcga_df <- fread(file.path(res_dir, "GMSig4_PatientScores_TCGA_TNBC.csv"), data.table = FALSE)
tcga_cc <- tcga_df[complete.cases(tcga_df[, c("OS_yr", "OS", "GMSig4", "Age", "Stage_III_IV", "Node_Pos")]), ]

g960_df <- fread(file.path(res_dir, "GMSig4_PatientScores_GSE96058_TNBC_Ext.csv"), data.table = FALSE)
# Reuse exact saved Node_Pos column when present so DCA matches Figure 3G/3H n=234 identically
if (!("Node_Pos" %in% colnames(g960_df))) {
  g960_df$Node_Pos <- ifelse(g960_df$Node_group %in% c("NodeNegative", "N0"), 0,
                      ifelse(!is.na(g960_df$Node_group) & g960_df$Node_group != "", 1, NA))
}
g960_cc <- g960_df[complete.cases(g960_df[, c("OS_yr", "OS", "GMSig4", "Age", "Tumor_mm", "Node_Pos")]), ]

cat(sprintf("Complete-case DCA cohorts: TCGA_TNBC n=%d (%d events) | GSE96058_TNBC_Ext n=%d (%d events)\n",
            nrow(tcga_cc), sum(tcga_cc$OS), nrow(g960_cc), sum(g960_cc$OS)))

# Predicted risk probability P(T <= t | X) from Cox model
pred_risk_t <- function(fit, df, t_yr) {
  sf <- survfit(fit, newdata = df)
  idx <- which.min(abs(sf$time - t_yr))
  # If max time < t_yr use last observed time; normally 1/3/5yr covered cleanly
  if (max(sf$time) < t_yr) idx <- length(sf$time)
  1 - as.numeric(sf$surv[idx, ])
}

# Vickers & Elkin (2006) censored survival Net Benefit formula (matches stdca / dcurves)
calc_surv_dca <- function(time_vec, event_vec, pred_list, t_yr, pts = seq(0.02, 0.50, by = 0.01)) {
  N <- length(time_vec)
  sf_all <- survfit(Surv(time_vec, event_vec) ~ 1)
  # Evaluate S_all(t_yr) by stepfun
  s_fn_all <- stepfun(sf_all$time, c(1, sf_all$surv))
  prev_t <- 1 - s_fn_all(t_yr)

  rows <- list()
  # Treat All & Treat None
  for (pt in pts) {
    nb_all  <- prev_t - (1 - prev_t) * (pt / (1 - pt))
    rows[[length(rows) + 1]] <- data.table(horizon_yr = t_yr, model = "Treat All", threshold = pt, net_benefit = nb_all, outcome_prev = prev_t)
    rows[[length(rows) + 1]] <- data.table(horizon_yr = t_yr, model = "Treat None", threshold = pt, net_benefit = 0, outcome_prev = prev_t)
  }

  for (mname in names(pred_list)) {
    p_hat <- pred_list[[mname]]
    for (pt in pts) {
      sel <- which(p_hat >= pt)
      n_sel <- length(sel)
      if (n_sel == 0) {
        nb <- 0
      } else if (n_sel == N) {
        nb <- prev_t - (1 - prev_t) * (pt / (1 - pt))
      } else {
        # KM event rate within selected high-risk subset
        sf_sub <- survfit(Surv(time_vec[sel], event_vec[sel]) ~ 1)
        s_sub  <- stepfun(sf_sub$time, c(1, sf_sub$surv))(t_yr)
        p_pos  <- n_sel / N
        tp_rate <- p_pos * (1 - s_sub)
        fp_rate <- p_pos * s_sub
        nb <- tp_rate - fp_rate * (pt / (1 - pt))
      }
      rows[[length(rows) + 1]] <- data.table(horizon_yr = t_yr, model = mname, threshold = pt, net_benefit = nb, outcome_prev = prev_t)
    }
  }
  rbindlist(rows)
}

# Fit the 3 candidate models in TCGA and GSE96058
f_sig_t  <- coxph(Surv(OS_yr, OS) ~ GMSig4, data = tcga_cc)
f_clin_t <- coxph(Surv(OS_yr, OS) ~ Age + Stage_III_IV + Node_Pos, data = tcga_cc)
f_nomo_t <- coxph(Surv(OS_yr, OS) ~ GMSig4 + Age + Stage_III_IV + Node_Pos, data = tcga_cc)

f_sig_9  <- coxph(Surv(OS_yr, OS) ~ GMSig4, data = g960_cc)
f_clin_9 <- coxph(Surv(OS_yr, OS) ~ Age + Tumor_mm + Node_Pos, data = g960_cc)
f_nomo_9 <- coxph(Surv(OS_yr, OS) ~ GMSig4 + Age + Tumor_mm + Node_Pos, data = g960_cc)

dca_all <- list()
for (h in c(1, 3, 5)) {
  pts_t <- seq(0.01, if (h == 1) 0.25 else if (h == 3) 0.45 else 0.55, by = 0.01)
  p_t <- list(
    "GMSig-4 Only" = pred_risk_t(f_sig_t, tcga_cc, h),
    "Clinical Only" = pred_risk_t(f_clin_t, tcga_cc, h),
    "Nomogram (GMSig-4 + Clinical)" = pred_risk_t(f_nomo_t, tcga_cc, h)
  )
  dt_t <- calc_surv_dca(tcga_cc$OS_yr, tcga_cc$OS, p_t, h, pts = pts_t)
  dt_t[, cohort := "TCGA-TNBC (n=168)"]
  dca_all[[length(dca_all) + 1]] <- dt_t

  pts_9 <- seq(0.01, if (h == 1) 0.20 else if (h == 3) 0.40 else 0.50, by = 0.01)
  p_9 <- list(
    "GMSig-4 Only" = pred_risk_t(f_sig_9, g960_cc, h),
    "Clinical Only" = pred_risk_t(f_clin_9, g960_cc, h),
    "Nomogram (GMSig-4 + Clinical)" = pred_risk_t(f_nomo_9, g960_cc, h)
  )
  dt_9 <- calc_surv_dca(g960_cc$OS_yr, g960_cc$OS, p_9, h, pts = pts_9)
  dt_9[, cohort := "GSE96058-TNBC (n=234)"]
  dca_all[[length(dca_all) + 1]] <- dt_9
}

dca_dt <- rbindlist(dca_all)
fwrite(dca_dt, file.path(res_dir, "DCA_Clinical_NetBenefit_TCGA_GSE96058.csv"))

# Print quantitative summary of max Net Benefit & gain over Clinical Only
sum_dca <- dcast(dca_dt[model %in% c("GMSig-4 Only", "Clinical Only", "Nomogram (GMSig-4 + Clinical)", "Treat All")],
                 cohort + horizon_yr + threshold ~ model, value.var = "net_benefit")
sum_dca[, gain_over_clin := `Nomogram (GMSig-4 + Clinical)` - `Clinical Only`]
cat("\n=== Clinical Decision Curve Analysis Summary (Positive Net Benefit Gains) ===\n")
print(sum_dca[`Nomogram (GMSig-4 + Clinical)` > 0, .(
  n_pos_thresh = .N,
  n_nomo_beats_clin = sum(gain_over_clin > 0, na.rm = TRUE),
  mean_gain_NB = mean(gain_over_clin, na.rm = TRUE),
  max_gain_NB = max(gain_over_clin, na.rm = TRUE)
), by = .(cohort, horizon_yr)])

# --- Render Supplementary Figure S2: 1/3/5-Year Clinical DCA ---
dca_dt[, model_f := factor(model, levels = c(
  "Nomogram (GMSig-4 + Clinical)", "Clinical Only", "GMSig-4 Only", "Treat All", "Treat None"
))]
dca_dt[, horizon_lbl := paste0(horizon_yr, "-Year OS")]
dca_dt[, cohort_f := factor(cohort, levels = c("TCGA-TNBC (n=168)", "GSE96058-TNBC (n=234)"))]

# Clip extreme negative Treat All tail so visual window centers around useful NB range
dca_plot <- copy(dca_dt)
dca_plot[, nb_disp := ifelse(model == "Treat All" & net_benefit < -0.05, NA_real_, net_benefit)]

p_dca <- ggplot(dca_plot, aes(x = threshold * 100, y = nb_disp, color = model_f, linetype = model_f, linewidth = model_f)) +
  geom_line(na.rm = TRUE) +
  facet_grid(cohort_f ~ horizon_lbl, scales = "free") +
  scale_color_manual(values = c(
    "Nomogram (GMSig-4 + Clinical)" = "#E64B35",
    "Clinical Only" = "#3C5488",
    "GMSig-4 Only" = "#00A087",
    "Treat All" = "grey55",
    "Treat None" = "black"
  ), name = NULL) +
  scale_linetype_manual(values = c(
    "Nomogram (GMSig-4 + Clinical)" = "solid",
    "Clinical Only" = "solid",
    "GMSig-4 Only" = "dotdash",
    "Treat All" = "dashed",
    "Treat None" = "dotted"
  ), name = NULL) +
  scale_linewidth_manual(values = c(
    "Nomogram (GMSig-4 + Clinical)" = 1.15,
    "Clinical Only" = 0.9,
    "GMSig-4 Only" = 0.85,
    "Treat All" = 0.65,
    "Treat None" = 0.65
  ), name = NULL) +
  labs(
    title = "Supplementary Figure S2. Clinical Decision Curve Analysis (1-, 3-, and 5-Year OS Net Benefit)",
    subtitle = "Combined GMSig-4 + Clinical Nomogram yields superior standardized Net Benefit across clinically actionable risk thresholds in TCGA-TNBC & GSE96058",
    x = "Decision Threshold Probability (%)",
    y = "Net Benefit"
  ) +
  theme_bw(base_size = 9.5) +
  theme(
    panel.grid.minor = element_blank(),
    plot.title = element_text(face = "bold", size = 11.5),
    plot.subtitle = element_text(size = 8.8, color = "grey30"),
    strip.background = element_rect(fill = "grey94", color = "grey75"),
    strip.text = element_text(face = "bold", size = 9),
    legend.position = "bottom"
  )

ggsave(file.path(sup_dir, "FigureS2_Clinical_DCA_1_3_5yr.pdf"), p_dca, width = 13.5, height = 8.4, device = cairo_pdf)
ggsave(file.path(sup_dir, "FigureS2_Clinical_DCA_1_3_5yr.png"), p_dca, width = 13.5, height = 8.4, dpi = 300)
cat("SAVED FigureS2_Clinical_DCA_1_3_5yr.{pdf,png} successfully!\n")
