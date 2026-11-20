project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(limma)
  library(survival)
})

data_dir <- file.path(project_root, "01_data", "TCGA_GEO")
res_dir  <- file.path(project_root, "02_results")
dir.create(res_dir, showWarnings = FALSE, recursive = TRUE)

cat("=== 1. Load TCGA-BRCA Clinical & Survival ===\n")
clin <- fread(file.path(data_dir, "TCGA_BRCA_clinicalMatrix.tsv"), data.table = FALSE)
surv <- fread(file.path(data_dir, "TCGA_BRCA_survival.tsv"), data.table = FALSE)

# Standardize TNBC definition via IHC / FISH
er  <- clin$breast_carcinoma_estrogen_receptor_status
pr  <- clin$breast_carcinoma_progesterone_receptor_status
her_ihc  <- clin$lab_proc_her2_neu_immunohistochemistry_receptor_status
her_fish <- clin$lab_procedure_her2_neu_in_situ_hybrid_outcome_type

# Fallback to nature2012 if blank
er  <- ifelse(is.na(er) | er == "", clin$ER_Status_nature2012, er)
pr  <- ifelse(is.na(pr) | pr == "", clin$PR_Status_nature2012, pr)
her <- ifelse(her_ihc %in% c("Positive", "Negative"), her_ihc,
              ifelse(her_fish %in% c("Positive", "Negative"), her_fish, clin$HER2_Final_Status_nature2012))

clin$ER_clean   <- er
clin$PR_clean   <- pr
clin$HER2_clean <- her
clin$is_TNBC    <- (clin$sample_type == "Primary Tumor") & (er == "Negative") & (pr == "Negative") & (her == "Negative")
clin$is_Normal  <- (clin$sample_type == "Solid Tissue Normal")
clin$is_NonTNBC <- (clin$sample_type == "Primary Tumor") & !clin$is_TNBC & (er %in% c("Positive", "Negative"))

cat(sprintf("Clinical counts -> Primary TNBC: %d | Primary Non-TNBC: %d | Normal: %d\n",
            sum(clin$is_TNBC, na.rm = TRUE), sum(clin$is_NonTNBC, na.rm = TRUE), sum(clin$is_Normal, na.rm = TRUE)))

cat("=== 2. Load TCGA-BRCA HiSeqV2 Expression ===\n")
expr_dt <- fread(file.path(data_dir, "TCGA_BRCA_HiSeqV2.tsv.gz"), data.table = FALSE)
genes <- expr_dt[[1]]
expr_mat <- as.matrix(expr_dt[, -1])
rownames(expr_mat) <- genes
rm(expr_dt); gc()

common_samples <- intersect(colnames(expr_mat), clin$sampleID)
clin <- clin[match(common_samples, clin$sampleID), ]
clin <- merge(clin, surv, by.x = "sampleID", by.y = "sample", all.x = TRUE, sort = FALSE)
clin <- clin[match(common_samples, clin$sampleID), ]

cat(sprintf("Expression matched samples -> Total: %d | TNBC: %d | Non-TNBC: %d | Normal: %d\n",
            length(common_samples), sum(clin$is_TNBC, na.rm=TRUE), sum(clin$is_NonTNBC, na.rm=TRUE), sum(clin$is_Normal, na.rm=TRUE)))

# 3. Differential expression: TNBC vs Normal & TNBC vs Non-TNBC
idx_tn_norm <- which(clin$is_TNBC | clin$is_Normal)
grp1 <- factor(ifelse(clin$is_TNBC[idx_tn_norm], "TNBC", "Normal"), levels = c("Normal", "TNBC"))
fit1 <- eBayes(contrasts.fit(lmFit(expr_mat[, idx_tn_norm], model.matrix(~ 0 + grp1)),
                             makeContrasts(grp1TNBC - grp1Normal, levels = model.matrix(~ 0 + grp1))))
deg_tn_norm <- topTable(fit1, number = Inf, adjust.method = "BH")
deg_tn_norm$gene <- rownames(deg_tn_norm)

idx_tn_non <- which(clin$is_TNBC | clin$is_NonTNBC)
grp2 <- factor(ifelse(clin$is_TNBC[idx_tn_non], "TNBC", "NonTNBC"), levels = c("NonTNBC", "TNBC"))
fit2 <- eBayes(contrasts.fit(lmFit(expr_mat[, idx_tn_non], model.matrix(~ 0 + grp2)),
                             makeContrasts(grp2TNBC - grp2NonTNBC, levels = model.matrix(~ 0 + grp2))))
deg_tn_non <- topTable(fit2, number = Inf, adjust.method = "BH")
deg_tn_non$gene <- rownames(deg_tn_non)

# 4. Intersect with GMrepo-anchored gutMGene targets & full gutMGene human targets
gut_gmrepo <- fread(file.path(res_dir, "GMrepo_anchored_gutMGene_targets.csv"), data.table = FALSE)
gut_master <- fread(file.path(res_dir, "gutMGene2_human_target_genes_master.csv"), data.table = FALSE)

# 5. Evaluate Univariate Cox (OS) & Immune Checkpoint / CD8 correlation within TNBC
idx_tnbc <- which(clin$is_TNBC)
tnbc_expr <- expr_mat[, idx_tnbc, drop = FALSE]
tnbc_clin <- clin[idx_tnbc, ]
valid_os <- which(!is.na(tnbc_clin$OS) & !is.na(tnbc_clin$OS.time) & tnbc_clin$OS.time > 0)
cat(sprintf("TNBC with valid OS > 0: %d (events=%d)\n", length(valid_os), sum(tnbc_clin$OS[valid_os] == 1)))

# Immune markers in TNBC
imm_markers <- c("CD8A", "GZMB", "PRF1", "IFNG", "CXCL13", "CD274", "PDCD1", "CTLA4")
imm_markers <- imm_markers[imm_markers %in% rownames(tnbc_expr)]
imm_score <- colMeans(tnbc_expr[imm_markers, , drop = FALSE])

target_genes <- intersect(gut_master$gene, rownames(expr_mat))
cat(sprintf("gutMGene human targets present in TCGA-BRCA: %d / %d (GMrepo-anchored: %d / %d)\n",
            length(target_genes), nrow(gut_master),
            sum(gut_gmrepo$gene %in% rownames(expr_mat)), nrow(gut_gmrepo)))

rows <- lapply(target_genes, function(g) {
  x <- as.numeric(tnbc_expr[g, ])
  # Cox OS
  cox_res <- tryCatch({
    fit <- coxph(Surv(OS.time, OS) ~ x_sub, data = data.frame(
      OS.time = tnbc_clin$OS.time[valid_os],
      OS = tnbc_clin$OS[valid_os],
      x_sub = x[valid_os]
    ))
    s <- summary(fit)
    c(HR = s$coefficients[1, "exp(coef)"], cox_p = s$coefficients[1, "Pr(>|z|)"])
  }, error = function(e) c(HR = NA, cox_p = NA))

  cor_imm <- cor.test(x, imm_score, method = "spearman")
  cor_pdl1 <- if ("CD274" %in% rownames(tnbc_expr)) cor.test(x, as.numeric(tnbc_expr["CD274", ]), method = "spearman") else list(estimate=NA, p.value=NA)
  cor_cd8 <- if ("CD8A" %in% rownames(tnbc_expr)) cor.test(x, as.numeric(tnbc_expr["CD8A", ]), method = "spearman") else list(estimate=NA, p.value=NA)

  r1 <- deg_tn_norm[g, ]
  r2 <- deg_tn_non[g, ]
  data.frame(
    gene = g,
    anchored_GMrepo = g %in% gut_gmrepo$gene,
    logFC_TNBC_vs_Normal = r1$logFC,
    FDR_TNBC_vs_Normal = r1$adj.P.Val,
    logFC_TNBC_vs_NonTNBC = r2$logFC,
    FDR_TNBC_vs_NonTNBC = r2$adj.P.Val,
    TNBC_OS_HR = cox_res["HR"],
    TNBC_OS_p = cox_res["cox_p"],
    cor_CytotoxicImmune_rho = as.numeric(cor_imm$estimate),
    cor_CytotoxicImmune_p = cor_imm$p.value,
    cor_CD8A_rho = as.numeric(cor_cd8$estimate),
    cor_CD274_rho = as.numeric(cor_pdl1$estimate),
    stringsAsFactors = FALSE
  )
})

df_screen <- rbindlist(rows)
df_screen <- merge(df_screen, gut_master[, c("gene", "n_microbes", "n_metabolites", "microbes", "metabolites")], by = "gene", all.x = TRUE)
df_screen <- df_screen[order(TNBC_OS_p, FDR_TNBC_vs_Normal)]
fwrite(df_screen, file.path(res_dir, "TCGA_TNBC_gutMGene_screen.csv"))

cat("\n=== Summary of Significant Hits (GMrepo-anchored & FDR<0.05 in TNBC vs Normal) ===\n")
deg_sig <- df_screen[anchored_GMrepo == TRUE & FDR_TNBC_vs_Normal < 0.05]
cat("GMrepo-anchored DEGs (TNBC vs Normal FDR<0.05):", nrow(deg_sig), "\n")
cat("  Also TNBC vs Non-TNBC FDR<0.05:", sum(deg_sig$FDR_TNBC_vs_NonTNBC < 0.05), "\n")
cat("  Prognostic in TNBC (OS p<0.05):", sum(df_screen$TNBC_OS_p < 0.05, na.rm=TRUE),
    "| GMrepo-anchored OS p<0.05:", sum(df_screen$anchored_GMrepo & df_screen$TNBC_OS_p < 0.05, na.rm=TRUE), "\n")
cat("  Strongly correlated with TNBC Cytotoxic/ICI score (|rho|>0.3, p<0.001):",
    sum(abs(df_screen$cor_CytotoxicImmune_rho) > 0.3 & df_screen$cor_CytotoxicImmune_p < 0.001, na.rm=TRUE), "\n\n")

print(head(df_screen[anchored_GMrepo == TRUE & FDR_TNBC_vs_Normal < 0.05 & (TNBC_OS_p < 0.1 | abs(cor_CytotoxicImmune_rho) > 0.3),
                     .(gene, logFC_TNBC_vs_Normal, FDR_TNBC_vs_Normal, logFC_TNBC_vs_NonTNBC, TNBC_OS_HR, TNBC_OS_p, cor_CytotoxicImmune_rho, cor_CD8A_rho, cor_CD274_rho)], 30))
