project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(survival)
  library(rms)
})

res_dir <- file.path(project_root, "02_results")
tcga_df <- fread(file.path(res_dir, "GMSig4_PatientScores_TCGA_TNBC.csv"), data.table = FALSE)
g960_df <- fread(file.path(res_dir, "GMSig4_PatientScores_GSE96058_TNBC_Ext.csv"), data.table = FALSE)

# Calibration points at 1, 3, 5 years using rms::cph + bootstrap calibrate (B=200)
compute_calib <- function(df, form_str, cohort_name, times_yr = c(1, 3, 5), m_grp = 45) {
  d <- df[complete.cases(df[, all.vars(as.formula(form_str))]), ]
  dd <- datadist(d)
  options(datadist = "dd")
  rows <- list()
  for (ty in times_yr) {
    if (sum(d$OS == 1 & d$OS_yr <= ty) < 5) next
    f <- cph(as.formula(form_str), data = d, x = TRUE, y = TRUE, surv = TRUE, time.inc = ty)
    cal <- tryCatch(calibrate(f, cmethod = "KM", method = "boot", u = ty, m = m_grp, B = 200), error = function(e) NULL)
    if (!is.null(cal)) {
      cal_df <- as.data.frame(cal[, c("mean.predicted", "KM", "std.err")])
      cal_df$cohort <- cohort_name
      cal_df$time_yr <- ty
      cal_df$n <- nrow(d)
      rows[[as.character(ty)]] <- cal_df
    }
  }
  rbindlist(rows, fill = TRUE)
}

# Set datadist globally cleanly inside environment for rms
dd_t <- datadist(tcga_df[complete.cases(tcga_df[, c("OS_yr", "OS", "GMSig4", "Age", "Stage_III_IV", "Node_Pos")]), ])
cal_t <- local({
  d <- tcga_df[complete.cases(tcga_df[, c("OS_yr", "OS", "GMSig4", "Age", "Stage_III_IV", "Node_Pos")]), ]
  dd <- datadist(d); assign("dd", dd, envir = .GlobalEnv); options(datadist = "dd")
  rows <- list()
  for (ty in c(1, 3, 5)) {
    f <- cph(Surv(OS_yr, OS) ~ GMSig4 + Age + Stage_III_IV + Node_Pos, data = d, x = TRUE, y = TRUE, surv = TRUE, time.inc = ty)
    cal <- calibrate(f, cmethod = "KM", method = "boot", u = ty, m = 42, B = 200)
    cd <- as.data.frame(cal[, c("mean.predicted", "KM", "std.err")])
    cd$cohort <- "TCGA_TNBC"; cd$time_yr <- ty; cd$n <- nrow(d)
    rows[[as.character(ty)]] <- cd
  }
  rbindlist(rows)
})

cal_9 <- local({
  d <- g960_df[complete.cases(g960_df[, c("OS_yr", "OS", "GMSig4", "Age", "Tumor_mm", "Node_Pos")]), ]
  dd <- datadist(d); assign("dd", dd, envir = .GlobalEnv); options(datadist = "dd")
  rows <- list()
  for (ty in c(1, 3, 5)) {
    f <- cph(Surv(OS_yr, OS) ~ GMSig4 + Age + Tumor_mm + Node_Pos, data = d, x = TRUE, y = TRUE, surv = TRUE, time.inc = ty)
    cal <- calibrate(f, cmethod = "KM", method = "boot", u = ty, m = 58, B = 200)
    cd <- as.data.frame(cal[, c("mean.predicted", "KM", "std.err")])
    cd$cohort <- "GSE96058_TNBC_Ext"; cd$time_yr <- ty; cd$n <- nrow(d)
    rows[[as.character(ty)]] <- cd
  }
  rbindlist(rows)
})

cal_all <- rbind(cal_t, cal_9)
fwrite(cal_all, file.path(res_dir, "GMSig4_Nomogram_Calibration_1_3_5yr.csv"))
print(cal_all)
cat("DONE Calibration computation!\n")
