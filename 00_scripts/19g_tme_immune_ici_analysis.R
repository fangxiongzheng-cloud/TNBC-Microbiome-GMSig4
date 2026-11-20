project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(GSVA)
  library(msigdbr)
  library(clusterProfiler)
  library(pROC)
})

cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")
res_dir   <- file.path(project_root, "02_results")

tcga <- readRDS(file.path(cache_dir, "TCGA_BRCA_clean.rds"))
g58  <- readRDS(file.path(cache_dir, "GSE58812_clean.rds"))
g960 <- readRDS(file.path(cache_dir, "GSE96058_TNBC_clean.rds"))
g194 <- readRDS(file.path(cache_dir, "GSE194040_clean.rds"))
g173 <- readRDS(file.path(cache_dir, "GSE173839_clean.rds"))

sc_tcga <- fread(file.path(res_dir, "GMSig4_PatientScores_TCGA_TNBC.csv"))
sc_g58  <- fread(file.path(res_dir, "GMSig4_PatientScores_GSE58812.csv"))
sc_g960 <- fread(file.path(res_dir, "GMSig4_PatientScores_GSE96058_TNBC_Ext.csv"))
mod_obj <- readRDS(file.path(res_dir, "GMSig4_Model_Objects.rds"))
w_ridge <- mod_obj$weights
w_uni   <- mod_obj$uni_betas; names(w_uni) <- names(w_ridge)

# 1. Curated 16 Immune / TME Functional Gene Sets + MSigDB Hallmarks
immune_sets <- list(
  CD8_Effector_T = c("CD8A", "CD8B", "GZMA", "GZMB", "PRF1", "IFNG", "TBX21", "EOMES", "CXCR3"),
  Cytolytic_Activity = c("GZMA", "PRF1", "GZMB", "GZMH", "GNLY"),
  Tcell_Inflamed_GEP = c("CCL5", "CD27", "CD274", "CD276", "CD8A", "CMKLR1", "CXCL9", "CXCR6", "HLA-DQA1", "HLA-DRB1", "HLA-E", "IDO1", "LAG3", "NKG7", "PDCD1LG2", "PSMB10", "STAT1", "TIGIT"),
  Tcell_Exhaustion = c("PDCD1", "CTLA4", "LAG3", "HAVCR2", "TIGIT", "TOX", "ENTPD1", "BTLA"),
  Antigen_Pres_MHCI = c("HLA-A", "HLA-B", "HLA-C", "B2M", "TAP1", "TAP2", "TAPBP", "PSMB8", "PSMB9"),
  Antigen_Pres_MHCII = c("HLA-DRA", "HLA-DRB1", "HLA-DPA1", "HLA-DPB1", "HLA-DQA1", "HLA-DQB1", "CIITA", "CD74"),
  M1_Macrophage = c("NOS2", "IL12A", "IL12B", "TNF", "IL1B", "CXCL9", "CXCL10", "CXCL11", "CD80", "CD86", "IRF5", "STAT1"),
  M2_Macrophage = c("CD163", "MRC1", "MSR1", "ARG1", "IL10", "TGFB1", "CCL18", "CCL22", "CSF1R"),
  Dendritic_Cells = c("BATF3", "IRF8", "CLEC9A", "XCR1", "CD1C", "FCER1A", "ITGAX", "LAMP3", "CCR7"),
  NK_Cells = c("NKG7", "GNLY", "KLRD1", "KLRK1", "KLRC1", "NCR1", "NCR3", "FCGR3A", "PRF1"),
  B_Cells = c("CD19", "MS4A1", "CD79A", "CD79B", "BLK", "BANK1", "PAX5", "CR2"),
  Treg = c("FOXP3", "IL2RA", "IKZF2", "TNFRSF18", "CTLA4", "CCR8", "TIGIT"),
  MDSC = c("S100A8", "S100A9", "S100A12", "ARG1", "IDO1", "NOS2", "CYBB", "CEACAM8", "ITGAM", "CD33"),
  CAF_ECM_Exclusion = c("FAP", "PDPN", "ACTA2", "COL1A1", "COL1A2", "COL3A1", "FN1", "POSTN", "MMP2", "TGFB1"),
  Glycolysis_Lactate = c("ENO1", "LDHA", "SLC2A1", "HK2", "PFKM", "PKM", "ALDOA", "PGK1", "GAPDH", "TPI1"),
  Microbial_SCFA_Sensing = c("FFAR2", "FFAR3", "HCAR2", "SLC5A8", "SLC16A1", "AHR", "NR1I2", "VDR", "TLR2", "TLR4", "NLRP6")
)

cat("=== 1. Loading Local MSigDB Hallmarks (H) & KEGG GMTs ===\n")
msig_dir <- file.path(project_root, "01_data", "TCGA_GEO", "msigdb")
h_df <- clusterProfiler::read.gmt(file.path(msig_dir, "hallmark.gmt"))
colnames(h_df) <- c("gs_name", "gene_symbol")
k_df <- clusterProfiler::read.gmt(file.path(msig_dir, "kegg.gmt"))
colnames(k_df) <- c("gs_name", "gene_symbol")
h_sets <- split(h_df$gene_symbol, h_df$gs_name)
all_sets <- c(immune_sets, h_sets)

run_ssgsea_safe <- function(mat, gsets) {
  gsets_f <- lapply(gsets, function(gs) intersect(gs, rownames(mat)))
  gsets_f <- gsets_f[sapply(gsets_f, length) >= 4]
  if (packageVersion("GSVA") >= "1.50.0") {
    p <- ssgseaParam(exprData = mat, geneSets = gsets_f, normalize = TRUE)
    gsva(p, verbose = FALSE)
  } else {
    gsva(mat, gsets_f, method = "ssgsea", ssgsea.norm = TRUE, verbose = FALSE)
  }
}

run_tme_cohort <- function(expr, sc_dt, coh_name) {
  mat <- expr[, sc_dt$sampleID, drop = FALSE]
  ss  <- run_ssgsea_safe(as.matrix(mat), all_sets)
  saveRDS(ss, file.path(res_dir, paste0("ssGSEA_", coh_name, ".rds")))
  # Compare High vs Low + Spearman vs GMSig4
  rows <- vector("list", nrow(ss))
  hi <- sc_dt$RiskGroup == "High"
  lo <- sc_dt$RiskGroup == "Low"
  for (i in seq_len(nrow(ss))) {
    v <- as.numeric(ss[i, ])
    ct <- cor.test(v, sc_dt$GMSig4, method = "spearman", exact = FALSE)
    wt <- wilcox.test(v[hi], v[lo])
    rows[[i]] <- data.table(
      cohort = coh_name, signature = rownames(ss)[i],
      category = ifelse(rownames(ss)[i] %in% names(immune_sets), "Immune_TME", "Hallmark"),
      mean_High = mean(v[hi]), mean_Low = mean(v[lo]),
      diff_High_minus_Low = mean(v[hi]) - mean(v[lo]),
      wilcox_p = wt$p.value,
      spearman_rho = as.numeric(ct$estimate), spearman_p = ct$p.value
    )
  }
  dt <- rbindlist(rows)
  dt[, fdr := p.adjust(spearman_p, method = "BH"), by = category]
  dt
}

cat("Running ssGSEA across TCGA-TNBC (n=172), GSE58812 (n=107), GSE96058-TNBC (n=242)...\n")
tme_t  <- run_tme_cohort(tcga$expr, sc_tcga, "TCGA_TNBC")
tme_58 <- run_tme_cohort(g58$expr,  sc_g58,  "GSE58812")
tme_96 <- run_tme_cohort(g960$expr, sc_g960, "GSE96058_TNBC_Ext")
tme_all <- rbind(tme_t, tme_58, tme_96)
fwrite(tme_all, file.path(res_dir, "GMSig4_TME_Hallmark_ssGSEA_3Cohorts.csv"))

# Print top replicated Immune_TME and Hallmark signatures
cat("\n=== Top Immune_TME Signatures Correlated with GMSig4 Across All 3 Cohorts ===\n")
imm_w <- dcast(tme_all[category == "Immune_TME"], signature ~ cohort, value.var = c("spearman_rho", "spearman_p"))
print(imm_w[order(spearman_rho_TCGA_TNBC)])

cat("\n=== Top Hallmark Pathways Correlated with GMSig4 Across All 3 Cohorts ===\n")
hal_w <- dcast(tme_all[category == "Hallmark"], signature ~ cohort, value.var = c("spearman_rho", "spearman_p"))
print(head(hal_w[order(abs(spearman_rho_TCGA_TNBC) + abs(spearman_rho_GSE58812) + abs(spearman_rho_GSE96058_TNBC_Ext), decreasing = TRUE)], 20))

# 2. Checkpoints & Cytotoxic Genes across 3 cohorts
icp_genes <- c("CD274", "PDCD1", "CTLA4", "LAG3", "TIGIT", "HAVCR2", "CD8A", "GZMB", "PRF1", "STAT1", "CCL5", "CXCL9", "HLA-A", "HLA-DRA")
icp_rows <- list()
for (coh in list(list("TCGA_TNBC", tcga$expr, sc_tcga),
                 list("GSE58812", g58$expr, sc_g58),
                 list("GSE96058_TNBC_Ext", g960$expr, sc_g960))) {
  cn <- coh[[1]]; ex <- coh[[2]]; scd <- coh[[3]]
  mat <- ex[, scd$sampleID, drop = FALSE]
  hi <- scd$RiskGroup == "High"; lo <- scd$RiskGroup == "Low"
  for (g in icp_genes) {
    if (!(g %in% rownames(mat))) next
    v <- as.numeric(mat[g, ])
    ct <- cor.test(v, scd$GMSig4, method = "spearman", exact = FALSE)
    wt <- wilcox.test(v[hi], v[lo])
    icp_rows[[length(icp_rows) + 1]] <- data.table(
      cohort = cn, gene = g, mean_High = mean(v[hi]), mean_Low = mean(v[lo]),
      wilcox_p = wt$p.value, spearman_rho = as.numeric(ct$estimate), spearman_p = ct$p.value
    )
  }
}
icp_dt <- rbindlist(icp_rows)
fwrite(icp_dt, file.path(res_dir, "GMSig4_ImmuneCheckpoints_3Cohorts.csv"))

# 3. Genome-wide limma DE + clusterProfiler GSEA on TCGA-TNBC High vs Low GMSig4
cat("\n=== 3. Genome-wide limma DE & GSEA (Hallmark + KEGG) in TCGA-TNBC ===\n")
mat_t <- as.matrix(tcga$expr[, sc_tcga$sampleID, drop = FALSE])
grp_f <- factor(sc_tcga$RiskGroup, levels = c("Low", "High"))
des   <- model.matrix(~ grp_f)
fit   <- eBayes(lmFit(mat_t, des))
tt    <- topTable(fit, coef = "grp_fHigh", number = Inf, sort.by = "t")
tt$gene <- rownames(tt)
fwrite(tt, file.path(res_dir, "GMSig4_TCGA_TNBC_HighVsLow_limma_DEGs.csv"))

gene_list <- sort(setNames(tt$t, tt$gene), decreasing = TRUE)
if (length(h_sets) > 0) {
  gsea_h <- GSEA(gene_list, TERM2GENE = h_df[, c("gs_name", "gene_symbol")], pvalueCutoff = 0.25, verbose = FALSE)
  fwrite(as.data.table(gsea_h@result), file.path(res_dir, "GMSig4_TCGA_TNBC_GSEA_Hallmarks.csv"))
  cat("Top GSEA Hallmarks (High vs Low GMSig4 in TCGA-TNBC):\n")
  print(head(as.data.table(gsea_h@result)[order(p.adjust), .(ID, NES, pvalue, p.adjust)], 15))
}
if (!is.null(k_df) && nrow(k_df) > 0) {
  gsea_k <- GSEA(gene_list, TERM2GENE = k_df[, c("gs_name", "gene_symbol")], pvalueCutoff = 0.25, verbose = FALSE)
  fwrite(as.data.table(gsea_k@result), file.path(res_dir, "GMSig4_TCGA_TNBC_GSEA_KEGG.csv"))
  cat("Top GSEA KEGG pathways (High vs Low GMSig4 in TCGA-TNBC):\n")
  print(head(as.data.table(gsea_k@result)[order(p.adjust), .(ID, NES, pvalue, p.adjust)], 15))
}

# 4. Comprehensive Immunotherapy pCR Validation + Specificity Interaction in GSE194040 & GSE173839
cat("\n=== 4. Immunotherapy Response Specificity & Interaction Analysis ===\n")
score_ext <- function(expr, idx, w = w_ridge) {
  Z <- scale(t(expr[names(w), idx, drop = FALSE]))
  as.numeric(Z %*% w)
}

p_all <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !is.na(g194$clin$pCR))
p_tn  <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
p_hrp <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !g194$clin$is_TNBC & !is.na(g194$clin$pCR))
c_tn  <- which(g194$clin$arm == "Paclitaxel" & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
c_all <- which(g194$clin$arm == "Paclitaxel" & !is.na(g194$clin$pCR))

g194_sub_list <- list(
  list("GSE194040_Pembro_TNBC", p_tn),
  list("GSE194040_Pembro_All", p_all),
  list("GSE194040_Pembro_HRpHER2n", p_hrp),
  list("GSE194040_Paclitaxel_Ctrl_TNBC", c_tn),
  list("GSE194040_Paclitaxel_Ctrl_All", c_all)
)

ici_rows <- list()
for (item in g194_sub_list) {
  nm <- item[[1]]; idx <- item[[2]]
  s_r <- score_ext(g194$expr, idx, w_ridge)
  s_u <- score_ext(g194$expr, idx, w_uni)
  y01 <- g194$clin$pCR[idx]
  rc_r <- roc(y01, -s_r, quiet = TRUE, direction = "<")
  rc_u <- roc(y01, -s_u, quiet = TRUE, direction = "<")
  ci_r <- as.numeric(ci.auc(rc_r))
  low_grp <- s_r < median(s_r)
  pcr_low  <- mean(y01[low_grp] == 1)
  pcr_high <- mean(y01[!low_grp] == 1)
  ft_p <- fisher.test(table(low_grp, y01))$p.value
  wp_r <- wilcox.test(s_r[y01 == 1], s_r[y01 == 0])$p.value
  wp_u <- wilcox.test(s_u[y01 == 1], s_u[y01 == 0])$p.value

  # Benchmark PD-L1 (CD274), PD-1 (PDCD1), STAT1
  cd274_auc <- as.numeric(auc(roc(y01, as.numeric(g194$expr["CD274", idx]), quiet = TRUE, direction = "<")))
  pdcd1_auc <- as.numeric(auc(roc(y01, as.numeric(g194$expr["PDCD1", idx]), quiet = TRUE, direction = "<")))
  stat1_auc <- as.numeric(auc(roc(y01, as.numeric(g194$expr["STAT1", idx]), quiet = TRUE, direction = "<")))

  ici_rows[[length(ici_rows) + 1]] <- data.table(
    cohort = nm, n = length(idx), pCR_n = sum(y01 == 1),
    GMSig4_AUC = ci_r[2], GMSig4_AUC_LCI = ci_r[1], GMSig4_AUC_UCI = ci_r[3],
    GMSig4_wilcox_p = wp_r, pCR_rate_LowRisk = pcr_low, pCR_rate_HighRisk = pcr_high, Fisher_p = ft_p,
    UniBeta_AUC = as.numeric(auc(rc_u)), UniBeta_wilcox_p = wp_u,
    CD274_AUC = cd274_auc, PDCD1_AUC = pdcd1_auc, STAT1_AUC = stat1_auc
  )
}

# Formal Interaction test (Pembrolizumab vs Paclitaxel control) in I-SPY2
idx_inter_tn <- c(p_tn, c_tn)
df_itn <- data.frame(
  pCR = g194$clin$pCR[idx_inter_tn],
  is_Pembro = as.numeric(grepl("Pembrolizumab", g194$clin$arm[idx_inter_tn], ignore.case = TRUE)),
  LowRisk_Protective = as.numeric(score_ext(g194$expr, idx_inter_tn, w_ridge) < median(score_ext(g194$expr, idx_inter_tn, w_ridge))),
  neg_GMSig4 = -score_ext(g194$expr, idx_inter_tn, w_ridge)
)
fit_itn_bin <- glm(pCR ~ LowRisk_Protective * is_Pembro, data = df_itn, family = binomial)
fit_itn_cont <- glm(pCR ~ neg_GMSig4 * is_Pembro, data = df_itn, family = binomial)
cat("I-SPY2 TNBC (Pembrolizumab vs Paclitaxel) Interaction GLM (LowRisk_Protective * is_Pembro):\n")
print(summary(fit_itn_bin)$coefficients)
cat("I-SPY2 TNBC Continuous Interaction GLM (neg_GMSig4 * is_Pembro):\n")
print(summary(fit_itn_cont)$coefficients)

idx_inter_all <- c(p_all, c_all)
df_iall <- data.frame(
  pCR = g194$clin$pCR[idx_inter_all],
  is_Pembro = as.numeric(grepl("Pembrolizumab", g194$clin$arm[idx_inter_all], ignore.case = TRUE)),
  LowRisk_Protective = as.numeric(score_ext(g194$expr, idx_inter_all, w_ridge) < median(score_ext(g194$expr, idx_inter_all, w_ridge))),
  neg_GMSig4 = -score_ext(g194$expr, idx_inter_all, w_ridge)
)
fit_iall_bin  <- glm(pCR ~ LowRisk_Protective * is_Pembro, data = df_iall, family = binomial)
fit_iall_cont <- glm(pCR ~ neg_GMSig4 * is_Pembro, data = df_iall, family = binomial)
cat("I-SPY2 ALL (Pembrolizumab vs Paclitaxel) Continuous Interaction GLM:\n")
print(summary(fit_iall_cont)$coefficients)

# GSE173839 Durvalumab+Olaparib
hr_col <- grep("^HR\\.status", colnames(g173$clin), value = TRUE)[1]
d_all <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1))
d_tn  <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0)
for (item in list(list("GSE173839_DurvaOlap_All", d_all), list("GSE173839_DurvaOlap_TNBC", d_tn))) {
  nm <- item[[1]]; idx <- item[[2]]
  s_r <- score_ext(g173$expr, idx, w_ridge)
  s_u <- score_ext(g173$expr, idx, w_uni)
  y01 <- g173$clin$pCR.status[idx]
  rc_r <- roc(y01, -s_r, quiet = TRUE, direction = "<")
  rc_u <- roc(y01, -s_u, quiet = TRUE, direction = "<")
  ci_r <- as.numeric(ci.auc(rc_r))
  low_grp <- s_r < median(s_r)
  ici_rows[[length(ici_rows) + 1]] <- data.table(
    cohort = nm, n = length(idx), pCR_n = sum(y01 == 1),
    GMSig4_AUC = ci_r[2], GMSig4_AUC_LCI = ci_r[1], GMSig4_AUC_UCI = ci_r[3],
    GMSig4_wilcox_p = wilcox.test(s_r[y01 == 1], s_r[y01 == 0])$p.value,
    pCR_rate_LowRisk = mean(y01[low_grp] == 1), pCR_rate_HighRisk = mean(y01[!low_grp] == 1),
    Fisher_p = fisher.test(table(low_grp, y01))$p.value,
    UniBeta_AUC = as.numeric(auc(rc_u)),
    UniBeta_wilcox_p = wilcox.test(s_u[y01 == 1], s_u[y01 == 0])$p.value,
    CD274_AUC = as.numeric(auc(roc(y01, as.numeric(g173$clin$PDL1[idx]), quiet = TRUE, direction = "<"))),
    PDCD1_AUC = as.numeric(auc(roc(y01, as.numeric(g173$clin$PD1[idx]), quiet = TRUE, direction = "<"))),
    STAT1_AUC = as.numeric(auc(roc(y01, as.numeric(g173$clin$STAT1_sig[idx]), quiet = TRUE, direction = "<")))
  )
}
ici_dt <- rbindlist(ici_rows)
print(ici_dt)
fwrite(ici_dt, file.path(res_dir, "GMSig4_Immunotherapy_pCR_Validation_Summary.csv"))

# Save individual ICI patient scores for violin/ROC plots in Figure 5
save_ici_pts <- function(expr, idx, y01, arm, tnbc_flag) {
  data.table(
    GMSig4 = score_ext(expr, idx, w_ridge),
    GMSig4_UniBeta = score_ext(expr, idx, w_uni),
    RiskGroup = ifelse(score_ext(expr, idx, w_ridge) > median(score_ext(expr, idx, w_ridge)), "High", "Low"),
    pCR = y01, arm = arm, is_TNBC = tnbc_flag,
    scale(t(expr[names(w_ridge), idx, drop = FALSE]))
  )
}
p_pts <- rbind(
  save_ici_pts(g194$expr, p_all, g194$clin$pCR[p_all], "Pembrolizumab", g194$clin$is_TNBC[p_all]),
  save_ici_pts(g194$expr, c_all, g194$clin$pCR[c_all], "Paclitaxel_Ctrl", g194$clin$is_TNBC[c_all])
)
fwrite(p_pts, file.path(res_dir, "GMSig4_PatientScores_GSE194040_ISPY2.csv"))

d_pts <- save_ici_pts(g173$expr, d_all, g173$clin$pCR.status[d_all], "Durvalumab_Olaparib", g173$clin[[hr_col]][d_all] == 0)
d_pts[, `:=`(PD1 = g173$clin$PD1[d_all], PDL1 = g173$clin$PDL1[d_all], TIS_sig = g173$clin$TIS_sig[d_all], STAT1_sig = g173$clin$STAT1_sig[d_all])]
fwrite(d_pts, file.path(res_dir, "GMSig4_PatientScores_GSE173839_Durva.csv"))
cat("DONE Step 4 TME + GSEA + Immunotherapy Specificity & pCR!\n")
