#!/usr/bin/env Rscript
# Stage B: Strictly Independent Locked Validation, Clinical Incremental Value,
# Head-to-Head Comparator Benchmarks & Ablation (Protocol v1.1).
# Uses frozen TCGA-TNBC Ridge coefficients from frozen_coefficients.json WITHOUT refitting.

.libPaths(c("/root/miniconda3/envs/bioinfo-r-core/lib/R/library", "/root/R/tnbc_microbiome_lib", .libPaths()))
suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(timeROC)
  library(pROC)
  library(jsonlite)
})

root     <- "/mnt/d/Bioinformatics/TNBC_Microbiome"
harm_dir <- file.path(root, "01_data/JTM_upgrade/external/harmonized")
out_dir  <- file.path(root, "02_results/JTM_upgrade/independent_validation")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

source(file.path(root, "00_scripts/26a_gmsig4_score.R"))
fz <- fromJSON(file.path(root, "02_results/JTM_upgrade/model_freeze/frozen_coefficients.json"))
w_ridge <- unlist(fz$weights)
w_imm3  <- w_ridge[c("IFNG", "IDO1", "CXCL10")]
gep18_genes <- c("CCL5","CD27","CD274","CD276","CD8A","CMKLR1","CXCL9","CXCR6",
                 "HLA-DQA1","HLA-DRB1","HLA-E","IDO1","LAG3","NKG7","PDCD1LG2",
                 "PSMB10","STAT1","TIGIT")

load_expr_mat <- function(path) {
  dt <- fread(path)
  genes <- dt[[1]]
  mat <- as.matrix(dt[, -1, with = FALSE])
  rownames(mat) <- genes
  storage.mode(mat) <- "numeric"
  mat
}

calc_c_index <- function(time_yr, event, score) {
  y <- Surv(time_yr, event)
  concordance(y ~ score, reverse = TRUE)$concordance
}

boot_c_ci <- function(time_yr, event, score, B = 1000L, seed = 20261122L) {
  set.seed(seed)
  n <- length(time_yr)
  est <- calc_c_index(time_yr, event, score)
  vals <- numeric(B)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    if (sum(event[ii]) < 2L || length(unique(score[ii])) < 2L) {
      vals[b] <- NA_real_
    } else {
      vals[b] <- tryCatch(calc_c_index(time_yr[ii], event[ii], score[ii]), error = function(e) NA_real_)
    }
  }
  ok <- vals[is.finite(vals)]
  list(est = est, lci = unname(quantile(ok, 0.025)), uci = unname(quantile(ok, 0.975)), n_ok = length(ok))
}

paired_boot_delta_c <- function(time_yr, event, s1, s2, B = 1000L, seed = 20261122L) {
  set.seed(seed)
  n <- length(time_yr)
  c1 <- calc_c_index(time_yr, event, s1)
  c2 <- calc_c_index(time_yr, event, s2)
  delta <- c1 - c2
  diffs <- numeric(B)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    if (sum(event[ii]) < 2L) {
      diffs[b] <- NA_real_
    } else {
      diffs[b] <- tryCatch(
        calc_c_index(time_yr[ii], event[ii], s1[ii]) - calc_c_index(time_yr[ii], event[ii], s2[ii]),
        error = function(e) NA_real_
      )
    }
  }
  ok <- diffs[is.finite(diffs)]
  p_two <- 2 * min(mean(ok <= 0), mean(ok >= 0))
  p_two <- min(1.0, max(p_two, 1 / (length(ok) + 1)))
  list(c1 = c1, c2 = c2, delta = delta,
       lci = unname(quantile(ok, 0.025)), uci = unname(quantile(ok, 0.975)),
       p_boot = p_two, n_ok = length(ok))
}

eval_endpoint <- function(label, tier, expr, clin, id_col, time_col, event_col, B = 1000L) {
  ok <- which(clin$has_all4_gmsig4 & !is.na(clin[[time_col]]) & !is.na(clin[[event_col]]) & clin[[time_col]] > 0)
  sub <- clin[ok]
  sids <- as.character(sub[[id_col]])
  sc_df <- score_gmsig4(expr, w_ridge, sids)
  raw_s <- sc_df$GMSig4
  z_s   <- as.numeric(scale(raw_s))
  ty    <- as.numeric(sub[[time_col]])
  ev    <- as.integer(sub[[event_col]])
  y     <- Surv(ty, ev)

  cx_raw <- coxph(y ~ raw_s)
  cx_sd  <- coxph(y ~ z_s)
  s_raw  <- summary(cx_raw)
  s_sd   <- summary(cx_sd)
  ph_p   <- tryCatch(cox.zph(cx_sd)$table[1, "p"], error = function(e) NA_real_)

  grp <- factor(ifelse(raw_s > median(raw_s), "High", "Low"), levels = c("Low", "High"))
  cx_g <- summary(coxph(y ~ grp))
  km_p <- 1 - pchisq(survdiff(y ~ grp)$chisq, 1)

  rkm <- survfit(Surv(ty, 1 - ev) ~ 1)
  rev_km <- unname(summary(rkm)$table["median"])
  ci_b <- boot_c_ci(ty, ev, raw_s, B = B)

  troc <- tryCatch(timeROC(T = ty, delta = ev, marker = raw_s, cause = 1,
                           weighting = "marginal", times = c(3, 5), iid = FALSE), error = function(e) NULL)
  auc3 <- if (!is.null(troc)) as.numeric(troc$AUC[which.min(abs(troc$times - 3))]) else NA_real_
  auc5 <- if (!is.null(troc)) as.numeric(troc$AUC[which.min(abs(troc$times - 5))]) else NA_real_

  # Build patient-level table with all comparators
  get_z <- function(g) if (g %in% rownames(expr) && all(is.finite(expr[g, sids]))) as.numeric(scale(expr[g, sids])) else rep(NA_real_, length(sids))
  z_ifng <- get_z("IFNG"); z_ido1 <- get_z("IDO1"); z_cx10 <- get_z("CXCL10"); z_eno1 <- get_z("ENO1")
  imm3_raw <- z_ifng * w_imm3["IFNG"] + z_ido1 * w_imm3["IDO1"] + z_cx10 * w_imm3["CXCL10"]
  cyt_raw  <- if (all(c("GZMA", "PRF1") %in% rownames(expr)) && all(is.finite(expr[c("GZMA", "PRF1"), sids]))) {
    -0.5 * (expr["GZMA", sids] + expr["PRF1", sids])
  } else rep(NA_real_, length(sids))

  gep_present <- intersect(gep18_genes, rownames(expr))
  gep_raw <- if (length(gep_present) >= 17L && all(is.finite(expr[gep_present, sids]))) {
    Z_gep <- scale(t(expr[gep_present, sids, drop = FALSE]))
    -rowMeans(Z_gep)
  } else rep(NA_real_, length(sids))

  pt_dt <- copy(sub)
  pt_dt[, `:=`(
    GMSig4 = raw_s, GMSig4_z = z_s, RiskGroup = as.character(grp),
    STAT1_risk = -get_z("STAT1"), PDCD1_risk = -get_z("PDCD1"), CD274_risk = -get_z("CD274"),
    CYT_risk = if (all(is.finite(cyt_raw))) as.numeric(scale(cyt_raw)) else NA_real_,
    GEP18_risk = if (all(is.finite(gep_raw))) as.numeric(scale(gep_raw)) else NA_real_,
    ENO1_only_risk = z_eno1,
    Immune3_frozen_risk = as.numeric(scale(imm3_raw))
  )]

  summary_row <- data.table(
    cohort_endpoint = label, tier = tier, n = length(ty), events = sum(ev),
    reverse_KM_followup_yr = rev_km, score_sd = sd(raw_s),
    HR_per_raw_unit = s_raw$conf.int[1, "exp(coef)"],
    HR_raw_LCI = s_raw$conf.int[1, "lower .95"], HR_raw_UCI = s_raw$conf.int[1, "upper .95"],
    HR_per_cohort_SD = s_sd$conf.int[1, "exp(coef)"],
    HR_SD_LCI = s_sd$conf.int[1, "lower .95"], HR_SD_UCI = s_sd$conf.int[1, "upper .95"],
    Cox_p = s_sd$coefficients[1, "Pr(>|z|)"], Schoenfeld_PH_p = ph_p,
    HighLow_HR = cx_g$conf.int[1, "exp(coef)"], HighLow_LCI = cx_g$conf.int[1, "lower .95"],
    HighLow_UCI = cx_g$conf.int[1, "upper .95"], KM_logrank_p = km_p,
    C_index = ci_b$est, C_LCI = ci_b$lci, C_UCI = ci_b$uci, boot_ok = ci_b$n_ok,
    AUC_3yr = auc3, AUC_5yr = auc5,
    n_at_risk_3yr = sum(ty >= 3), n_at_risk_5yr = sum(ty >= 5)
  )

  # Comparator & ablation benchmarks
  comp_cols <- c("STAT1_risk", "PDCD1_risk", "CD274_risk", "CYT_risk", "GEP18_risk", "ENO1_only_risk", "Immune3_frozen_risk")
  comp_names <- c("STAT1_single", "PDCD1_single", "CD274_single", "CYT_GZMA_PRF1", "GEP18_unweighted", "ENO1_only", "Immune3_frozen")
  comp_rows <- list()
  for (k in seq_along(comp_cols)) {
    v <- pt_dt[[comp_cols[k]]]
    if (all(is.na(v))) next
    cc <- which(is.finite(v) & is.finite(raw_s))
    cx_c <- summary(coxph(Surv(ty[cc], ev[cc]) ~ v[cc]))
    pb <- paired_boot_delta_c(ty[cc], ev[cc], raw_s[cc], v[cc], B = B)
    comp_rows[[length(comp_rows) + 1L]] <- data.table(
      cohort_endpoint = label, comparator = comp_names[k], n = length(cc), events = sum(ev[cc]),
      comp_HR_per_SD = cx_c$conf.int[1, "exp(coef)"],
      comp_HR_LCI = cx_c$conf.int[1, "lower .95"], comp_HR_UCI = cx_c$conf.int[1, "upper .95"],
      comp_Cox_p = cx_c$coefficients[1, "Pr(>|z|)"],
      GMSig4_C = pb$c1, Comp_C = pb$c2,
      Delta_C_GMSig4_minus_Comp = pb$delta,
      Delta_C_LCI = pb$lci, Delta_C_UCI = pb$uci,
      p_boot_two_sided = pb$p_boot, boot_ok = pb$n_ok
    )
  }
  comp_dt <- rbindlist(comp_rows)
  comp_dt[, p_holm := p.adjust(p_boot_two_sided, method = "holm")]

  list(summary = summary_row, comparators = comp_dt, patient_df = pt_dt)
}

# Optimism-corrected C-index via bootstrap refitting for cohort-fitted clinical+score models
boot_optimism_c <- function(df, time_col, event_col, rhs_formula, B = 500L, seed = 20261122L) {
  set.seed(seed)
  fml <- as.formula(paste0("Surv(", time_col, ", ", event_col, ") ~ ", rhs_formula))
  fit0 <- coxph(fml, data = df)
  lp0  <- predict(fit0, type = "lp")
  c_app <- calc_c_index(df[[time_col]], df[[event_col]], lp0)
  n <- nrow(df)
  opts <- numeric(B)
  for (b in seq_len(B)) {
    ii <- sample.int(n, n, replace = TRUE)
    dfb <- df[ii]
    fb <- tryCatch(coxph(fml, data = dfb), error = function(e) NULL)
    if (is.null(fb) || any(!is.finite(coef(fb)))) { opts[b] <- NA_real_; next }
    c_boot <- tryCatch(calc_c_index(dfb[[time_col]], dfb[[event_col]], predict(fb, newdata = dfb, type = "lp")), error = function(e) NA_real_)
    c_orig <- tryCatch(calc_c_index(df[[time_col]],  df[[event_col]],  predict(fb, newdata = df,  type = "lp")), error = function(e) NA_real_)
    opts[b] <- c_boot - c_orig
  }
  ok <- opts[is.finite(opts)]
  list(c_app = c_app, optimism = mean(ok), c_corr = c_app - mean(ok), n_ok = length(ok), fit = fit0)
}

# --- Load harmonized cohorts ---
mb_clin <- fread(file.path(harm_dir, "METABRIC_TNBC_clinical.csv"))
mb_expr <- load_expr_mat(file.path(harm_dir, "METABRIC_TNBC_expr.csv"))

g25_clin <- fread(file.path(harm_dir, "GSE25066_TNBC_clinical.csv"))
g25_med  <- load_expr_mat(file.path(harm_dir, "GSE25066_TNBC_expr_median.csv"))
g25_mean <- load_expr_mat(file.path(harm_dir, "GSE25066_TNBC_expr_mean.csv"))

g21_clin <- fread(file.path(harm_dir, "GSE21653_TNBC_clinical.csv"))
g21_med  <- load_expr_mat(file.path(harm_dir, "GSE21653_TNBC_expr_median.csv"))
g21_mean <- load_expr_mat(file.path(harm_dir, "GSE21653_TNBC_expr_mean.csv"))

cat("=== Running Stage B Strictly Independent Locked Validation ===\n")
res_mb_os   <- eval_endpoint("METABRIC_TNBC_OS",        "PRIMARY_STRICT_INDEP",   mb_expr,  mb_clin,  "SAMPLE_ID", "OS_yr",   "OS_event")
res_mb_erihc<- eval_endpoint("METABRIC_TNBC_OS_ER_IHC", "SENSITIVITY_ER_IHC",     mb_expr,  mb_clin[is_TNBC_erihc == TRUE], "SAMPLE_ID", "OS_yr", "OS_event")
res_mb_dss  <- eval_endpoint("METABRIC_TNBC_DSS",       "SECONDARY_STRICT_INDEP", mb_expr,  mb_clin,  "SAMPLE_ID", "OS_yr",   "DSS_event")
res_mb_rfs  <- eval_endpoint("METABRIC_TNBC_RFS",       "SECONDARY_STRICT_INDEP", mb_expr,  mb_clin,  "SAMPLE_ID", "RFS_yr",  "RFS_event")
res_25_med  <- eval_endpoint("GSE25066_TNBC_DRFS",      "SECONDARY_STRICT_INDEP", g25_med,  g25_clin, "sampleID",  "DRFS_yr", "DRFS_event")
res_25_mean <- eval_endpoint("GSE25066_TNBC_DRFS_mean", "SENSITIVITY_MEAN_PROBE", g25_mean, g25_clin, "sampleID",  "DRFS_yr", "DRFS_event")
res_21_med  <- eval_endpoint("GSE21653_TNBC_DFS",       "SECONDARY_STRICT_INDEP", g21_med,  g21_clin, "sampleID",  "DFS_yr",  "DFS_event")
res_21_mean <- eval_endpoint("GSE21653_TNBC_DFS_mean",  "SENSITIVITY_MEAN_PROBE", g21_mean, g21_clin, "sampleID",  "DFS_yr",  "DFS_event")

all_sum <- rbindlist(list(
  res_mb_os$summary, res_mb_erihc$summary, res_mb_dss$summary, res_mb_rfs$summary,
  res_25_med$summary, res_25_mean$summary, res_21_med$summary, res_21_mean$summary
))
all_comp <- rbindlist(list(
  res_mb_os$comparators, res_mb_dss$comparators, res_mb_rfs$comparators,
  res_25_med$comparators, res_21_med$comparators
))

fwrite(all_sum,  file.path(out_dir, "Independent_Validation_Summary.csv"))
fwrite(all_comp, file.path(out_dir, "Independent_Comparator_Benchmarks.csv"))
fwrite(res_mb_os$patient_df,  file.path(out_dir, "PatientScores_METABRIC_TNBC_OS.csv"))
fwrite(res_25_med$patient_df, file.path(out_dir, "PatientScores_GSE25066_TNBC_DRFS.csv"))
fwrite(res_21_med$patient_df, file.path(out_dir, "PatientScores_GSE21653_TNBC_DFS.csv"))

print(all_sum[, .(cohort_endpoint, n, events, HR_SD = sprintf("%.3f (%.2f-%.2f)", HR_per_cohort_SD, HR_SD_LCI, HR_SD_UCI),
                  Cox_p = signif(Cox_p, 3), KM_p = signif(KM_logrank_p, 3),
                  C = sprintf("%.3f (%.3f-%.3f)", C_index, C_LCI, C_UCI), AUC3 = round(AUC_3yr, 3), AUC5 = round(AUC_5yr, 3))])

cat("\n=== Comparator Benchmarks (METABRIC OS, GSE25066 DRFS, GSE21653 DFS) ===\n")
print(all_comp[cohort_endpoint %in% c("METABRIC_TNBC_OS", "GSE25066_TNBC_DRFS", "GSE21653_TNBC_DFS"),
               .(cohort_endpoint, comparator, comp_HR = round(comp_HR_per_SD, 3), comp_p = signif(comp_Cox_p, 3),
                 GMSig4_C = round(GMSig4_C, 3), Comp_C = round(Comp_C, 3),
                 Delta_C = sprintf("%+.3f (%.3f to %.3f)", Delta_C_GMSig4_minus_Comp, Delta_C_LCI, Delta_C_UCI),
                 p_boot = signif(p_boot_two_sided, 3), p_holm = signif(p_holm, 3))])

# --- Multivariable clinical adjustment & optimism-corrected cohort-fitted incremental value ---
cat("\n=== Multivariable Clinical Adjustment & Cohort-Fitted Incremental Value ===\n")
mv_rows <- list()
inc_rows <- list()

run_mv_block <- function(df, coh_lbl, time_col, event_col, covars, model_note) {
  d <- copy(df[complete.cases(df[, c(time_col, event_col, "GMSig4_z", covars), with = FALSE])])
  if ("Age" %in% covars) d[, Age10 := Age / 10]
  if ("Tumor_mm" %in% covars) d[, Size10 := Tumor_mm / 10]
  rhs_c <- paste(sapply(covars, function(x) if (x == "Age") "Age10" else if (x == "Tumor_mm") "Size10" else x), collapse = " + ")
  rhs_f <- paste("GMSig4_z +", rhs_c)
  b_c <- boot_optimism_c(d, time_col, event_col, rhs_c, B = 500L)
  b_f <- boot_optimism_c(d, time_col, event_col, rhs_f, B = 500L)
  lrt_p <- anova(b_c$fit, b_f$fit)[2, "Pr(>|Chi|)"]
  s <- summary(b_f$fit)
  epv <- sum(d[[event_col]]) / length(coef(b_f$fit))
  for (v in rownames(s$coefficients)) {
    mv_rows[[length(mv_rows) + 1L]] <<- data.table(
      cohort_model = coh_lbl, note = model_note, n = nrow(d), events = sum(d[[event_col]]),
      EPV = round(epv, 2), variable = v,
      HR = s$conf.int[v, "exp(coef)"], LCI = s$conf.int[v, "lower .95"],
      UCI = s$conf.int[v, "upper .95"], pval = s$coefficients[v, "Pr(>|z|)"]
    )
  }
  inc_rows[[length(inc_rows) + 1L]] <<- data.table(
    cohort_model = coh_lbl, note = model_note, n = nrow(d), events = sum(d[[event_col]]), EPV = round(epv, 2),
    C_clin_apparent = b_c$c_app, C_clin_opt_corr = b_c$c_corr,
    C_full_apparent = b_f$c_app, C_full_opt_corr = b_f$c_corr,
    Delta_C_opt_corr = b_f$c_corr - b_c$c_corr, LRT_p = lrt_p
  )
}

run_mv_block(res_mb_os$patient_df,  "METABRIC_TNBC_OS_PrimaryAdj",   "OS_yr",   "OS_event",   c("Age", "Tumor_mm", "Node_Pos"), "Primary clinical adjustment (Age/10yr + Size/10mm + Node+)")
run_mv_block(res_mb_os$patient_df,  "METABRIC_TNBC_OS_StageSens",    "OS_yr",   "OS_event",   c("Age", "Stage_III_IV"),         "Stage sensitivity (Age/10yr + Stage III/IV)")
run_mv_block(res_mb_dss$patient_df, "METABRIC_TNBC_DSS_PrimaryAdj",  "OS_yr",   "DSS_event",  c("Age", "Tumor_mm", "Node_Pos"), "Secondary DSS adjustment")
run_mv_block(res_mb_rfs$patient_df, "METABRIC_TNBC_RFS_PrimaryAdj",  "RFS_yr",  "RFS_event",  c("Age", "Tumor_mm", "Node_Pos"), "Secondary RFS adjustment")
run_mv_block(res_25_med$patient_df, "GSE25066_TNBC_DRFS_PrimaryAdj", "DRFS_yr", "DRFS_event", c("Age", "Node_Pos", "Stage_III_IV"), "Secondary DRFS adjustment")
run_mv_block(res_21_med$patient_df, "GSE21653_TNBC_DFS_Exploratory", "DFS_yr",  "DFS_event",  c("Age", "Node_Pos"),             "Exploratory (EPV < 10)")

# Also test incremental value over Immune3 and GEP18 in METABRIC OS!
mb_p <- res_mb_os$patient_df[complete.cases(res_mb_os$patient_df[, .(OS_yr, OS_event, GMSig4_z, Immune3_frozen_risk, GEP18_risk, Age, Tumor_mm, Node_Pos)])]
mb_p[, `:=`(Age10 = Age / 10, Size10 = Tumor_mm / 10)]
f_cl_im3 <- coxph(Surv(OS_yr, OS_event) ~ Immune3_frozen_risk + Age10 + Size10 + Node_Pos, data = mb_p)
f_cl_im3_eno1 <- coxph(Surv(OS_yr, OS_event) ~ Immune3_frozen_risk + ENO1_only_risk + Age10 + Size10 + Node_Pos, data = mb_p)
f_cl_gep <- coxph(Surv(OS_yr, OS_event) ~ GEP18_risk + Age10 + Size10 + Node_Pos, data = mb_p)
f_cl_gep_gm4 <- coxph(Surv(OS_yr, OS_event) ~ GEP18_risk + GMSig4_z + Age10 + Size10 + Node_Pos, data = mb_p)
cat("METABRIC OS LRT adding ENO1 to Clinical+Immune3: p =", signif(anova(f_cl_im3, f_cl_im3_eno1)[2, "Pr(>|Chi|)"], 4), "\n")
cat("METABRIC OS LRT adding GMSig4 to Clinical+GEP18: p =", signif(anova(f_cl_gep, f_cl_gep_gm4)[2, "Pr(>|Chi|)"], 4), "\n")

mv_dt  <- rbindlist(mv_rows)
inc_dt <- rbindlist(inc_rows)
fwrite(mv_dt,  file.path(out_dir, "Independent_Multivariable_Cox.csv"))
fwrite(inc_dt, file.path(out_dir, "Independent_Clinical_Incremental_Value.csv"))
print(inc_dt[, .(cohort_model, n, events, EPV, C_clin_opt = round(C_clin_opt_corr, 3), C_full_opt = round(C_full_opt_corr, 3), Delta_opt = round(Delta_C_opt_corr, 4), LRT_p = signif(LRT_p, 3))])

# Exploratory chemo pCR specificity check in GSE25066 (neoadjuvant chemo-only!)
pcr_sub <- res_25_med$patient_df[!is.na(pCR)]
rc_pcr <- roc(pcr_sub$pCR, -pcr_sub$GMSig4, quiet = TRUE, direction = "<")
ci_pcr <- as.numeric(ci.auc(rc_pcr))
cat(sprintf("\nExploratory GSE25066 Chemo-Only pCR (n=%d, pCR=%d): GMSig4 AUC = %.3f (95%% CI %.3f-%.3f), Wilcoxon p = %.4g\n",
            nrow(pcr_sub), sum(pcr_sub$pCR == 1), ci_pcr[2], ci_pcr[1], ci_pcr[3],
            wilcox.test(pcr_sub$GMSig4[pcr_sub$pCR == 1], pcr_sub$GMSig4[pcr_sub$pCR == 0])$p.value))
cat("STAGE B COMPLETE.\n")
