project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ggrepel)
  library(fgsea)
  library(clusterProfiler)
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

res_dir  <- file.path(project_root, "02_results")
fig_dir  <- file.path(project_root, "03_figures")
msig_dir <- file.path(project_root, "01_data", "TCGA_GEO", "msigdb")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

PAL <- list(
  low    = "#2B7A9E",
  high   = "#C84B31",
  emer   = "#1F8A70",
  navy   = "#1E3A5F",
  ochre  = "#D9822B",
  div_lo = "#1D5B9B",
  div_mi = "#F8FAFC",
  div_hi = "#C0392B",
  ink    = "#1A202C",
  muted  = "#4A5568",
  line   = "#2D3748"
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
      plot.margin = margin(6, 6, 4, 6)
    )
}

# ==============================================================================
# PANEL A: Full 18-Signature TME & Hallmark Correlation Heatmap w/ Module Strip
# ==============================================================================
ss_cor <- fread(file.path(res_dir, "GMSig4_TME_Hallmark_ssGSEA_3Cohorts.csv"))
key_sigs <- c(
  "Tcell_Inflamed_GEP",
  "CD8_Effector_T",
  "Cytolytic_Activity",
  "Antigen_Pres_MHCI",
  "Antigen_Pres_MHCII",
  "NK_Cells",
  "M1_Macrophage",
  "Dendritic_Cells",
  "Tcell_Exhaustion",
  "Microbial_SCFA_Sensing",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE",
  "HALLMARK_INFLAMMATORY_RESPONSE",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING",
  "CAF_ECM_Exclusion",
  "HALLMARK_ANGIOGENESIS",
  "HALLMARK_HYPOXIA",
  "HALLMARK_GLYCOLYSIS",
  "Glycolysis_Lactate",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION"
)
sig_labels <- c(
  "Tcell_Inflamed_GEP"                         = "T-cell Inflamed GEP",
  "CD8_Effector_T"                             = "CD8+ Effector T Cells",
  "Cytolytic_Activity"                         = "Cytolytic Activity",
  "Antigen_Pres_MHCI"                          = "Antigen Pres. (MHC-I)",
  "Antigen_Pres_MHCII"                         = "Antigen Pres. (MHC-II)",
  "NK_Cells"                                   = "NK Cell Cytotoxicity",
  "M1_Macrophage"                              = "M1 Macrophage",
  "Dendritic_Cells"                            = "Dendritic Cells",
  "Tcell_Exhaustion"                           = "T-cell Exhaustion",
  "Microbial_SCFA_Sensing"                     = "Microbial SCFA Sensing",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"         = "Hallmark: IFN-\u03B3 Response",
  "HALLMARK_INFLAMMATORY_RESPONSE"             = "Hallmark: Inflammatory Resp.",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"           = "Hallmark: IL6-JAK-STAT3",
  "CAF_ECM_Exclusion"                          = "CAF ECM Exclusion",
  "HALLMARK_ANGIOGENESIS"                      = "Hallmark: Angiogenesis",
  "HALLMARK_HYPOXIA"                           = "Hallmark: Hypoxia",
  "HALLMARK_GLYCOLYSIS"                        = "Hallmark: Glycolysis",
  "Glycolysis_Lactate"                         = "Glycolysis\u2013Lactate Axis",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Hallmark: EMT"
)
sig_mod <- c(
  "Tcell_Inflamed_GEP"                         = "Cytotoxic TME",
  "CD8_Effector_T"                             = "Cytotoxic TME",
  "Cytolytic_Activity"                         = "Cytotoxic TME",
  "Antigen_Pres_MHCI"                          = "Cytotoxic TME",
  "Antigen_Pres_MHCII"                         = "Cytotoxic TME",
  "NK_Cells"                                   = "Cytotoxic TME",
  "M1_Macrophage"                              = "Cytotoxic TME",
  "Dendritic_Cells"                            = "Cytotoxic TME",
  "Tcell_Exhaustion"                           = "Cytotoxic TME",
  "Microbial_SCFA_Sensing"                     = "Microbial / IFN",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"         = "Microbial / IFN",
  "HALLMARK_INFLAMMATORY_RESPONSE"             = "Microbial / IFN",
  "HALLMARK_IL6_JAK_STAT3_SIGNALING"           = "Microbial / IFN",
  "CAF_ECM_Exclusion"                          = "Stromal / Metabolic",
  "HALLMARK_ANGIOGENESIS"                      = "Stromal / Metabolic",
  "HALLMARK_HYPOXIA"                           = "Stromal / Metabolic",
  "HALLMARK_GLYCOLYSIS"                        = "Stromal / Metabolic",
  "Glycolysis_Lactate"                         = "Stromal / Metabolic",
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "Stromal / Metabolic"
)
mod_cols <- c("Cytotoxic TME" = PAL$low, "Microbial / IFN" = PAL$emer, "Stromal / Metabolic" = PAL$high)

sub_cor <- ss_cor[signature %in% key_sigs]
sub_cor[, sig_clean := factor(sig_labels[signature], levels = rev(unname(sig_labels[key_sigs])))]
sub_cor[, cohort_lbl := factor(cohort,
  levels = c("TCGA_TNBC", "GSE58812", "GSE96058_TNBC_Ext"),
  labels = c("TCGA-TNBC\n(n=172)", "GSE58812\n(n=107)", "GSE96058\n(n=242)")
)]
sub_cor[, stars := fifelse(spearman_p < 0.001, "***", fifelse(spearman_p < 0.01, "**", fifelse(spearman_p < 0.05, "*", "")))]
sub_cor[, txt_col := ifelse(abs(spearman_rho) > 0.50, "white", PAL$ink)]

strip_4a <- data.table(
  sig_clean = factor(rev(unname(sig_labels[key_sigs])), levels = rev(unname(sig_labels[key_sigs]))),
  Module = factor(sig_mod[rev(key_sigs)], levels = names(mod_cols))
)

p4a_strip <- ggplot(strip_4a, aes(x = 1, y = sig_clean, fill = Module)) +
  geom_tile(color = "white", linewidth = 0.45, width = 0.85) +
  scale_fill_manual(values = mod_cols, name = NULL) +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7.2, face = "bold"),
    legend.position = "bottom", legend.key.size = unit(6.5, "pt"),
    plot.margin = margin(6, 0, 4, 4)
  ) +
  guides(fill = guide_legend(nrow = 2))

p4a_mat <- ggplot(sub_cor, aes(x = cohort_lbl, y = sig_clean, fill = spearman_rho)) +
  geom_tile(color = "white", linewidth = 0.55) +
  geom_text(aes(label = sprintf("%.2f%s", spearman_rho, stars), color = txt_col),
            size = 2.35, fontface = "bold", family = "Arial") +
  scale_color_identity() +
  scale_fill_gradient2(low = PAL$div_lo, mid = PAL$div_mi, high = PAL$div_hi, midpoint = 0,
                       limits = c(-0.92, 0.55), name = "Spearman \u03C1") +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_blank(), axis.text.x = element_text(face = "bold", size = 7.2),
    legend.position = "bottom", plot.margin = margin(6, 6, 4, 1)
  ) +
  guides(fill = guide_colorbar(barwidth = unit(72, "pt"), barheight = unit(5.5, "pt"),
                               title.position = "left", title.vjust = 0.85))

p4a <- (p4a_strip | p4a_mat) + plot_layout(widths = c(0.09, 1.0), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(6, "pt"))

# ==============================================================================
# PANEL B: High- vs Low-GMSig4 TME Rainclouds (Half-Violin + Box + Dots)
# ==============================================================================
load_ssgsea_long <- function(rds_name, csv_name, coh_label) {
  mat <- as.matrix(readRDS(file.path(res_dir, rds_name)))
  sc  <- fread(file.path(res_dir, csv_name))
  axes <- c("CD8_Effector_T", "Tcell_Inflamed_GEP", "Microbial_SCFA_Sensing", "Glycolysis_Lactate")
  axes <- intersect(axes, rownames(mat))
  common <- intersect(sc$sampleID, colnames(mat))
  sc <- sc[match(common, sampleID)]
  dt_m <- as.data.table(t(mat[axes, common, drop = FALSE]))
  dt_m[, RiskGroup := sc$RiskGroup]
  dt_m[, Cohort := coh_label]
  for (ax in axes) dt_m[[ax]] <- as.numeric(scale(dt_m[[ax]]))
  melt(dt_m, id.vars = c("Cohort", "RiskGroup"), variable.name = "Axis", value.name = "Zscore")
}

b_dt <- rbindlist(list(
  load_ssgsea_long("ssGSEA_TCGA_TNBC.rds", "GMSig4_PatientScores_TCGA_TNBC.csv", "TCGA-TNBC"),
  load_ssgsea_long("ssGSEA_GSE58812.rds", "GMSig4_PatientScores_GSE58812.csv", "GSE58812"),
  load_ssgsea_long("ssGSEA_GSE96058_TNBC_Ext.rds", "GMSig4_PatientScores_GSE96058_TNBC_Ext.csv", "GSE96058")
))
axis_map <- c(
  "CD8_Effector_T"         = "CD8+ Effector T",
  "Tcell_Inflamed_GEP"     = "T-cell Inflamed GEP",
  "Microbial_SCFA_Sensing" = "Microbial SCFA Sensing",
  "Glycolysis_Lactate"     = "Glycolysis\u2013Lactate"
)
b_dt[, Axis_lbl := factor(axis_map[as.character(Axis)], levels = unname(axis_map))]
b_dt[, RiskGroup := factor(RiskGroup, levels = c("Low", "High"))]
b_dt[, Cohort := factor(Cohort, levels = c("TCGA-TNBC", "GSE58812", "GSE96058"))]

pvals_4b <- b_dt[, .(
  wilcox_p = wilcox.test(Zscore ~ RiskGroup)$p.value,
  y_pos = max(Zscore, na.rm = TRUE) + 0.40
), by = .(Cohort, Axis_lbl)]
pvals_4b[, p_lbl := ifelse(wilcox_p < 1e-4, sprintf("P=%.1e", wilcox_p), sprintf("P=%.4f", wilcox_p))]

p4b <- ggplot(b_dt, aes(x = Cohort, y = Zscore, fill = RiskGroup, color = RiskGroup)) +
  geom_violin(position = position_dodge(0.78), width = 0.72, alpha = 0.25, linewidth = 0.35, trim = TRUE) +
  geom_boxplot(position = position_dodge(0.78), width = 0.24, outlier.size = 0.45, outlier.alpha = 0.45,
               alpha = 0.85, color = PAL$ink, linewidth = 0.35) +
  geom_text(data = pvals_4b, aes(x = Cohort, y = y_pos, label = p_lbl),
            inherit.aes = FALSE, size = 2.15, fontface = "bold", family = "Arial", color = PAL$ink) +
  facet_wrap(~Axis_lbl, scales = "free_y", ncol = 2) +
  scale_fill_manual(values = c("Low" = PAL$low, "High" = PAL$high), name = "GMSig-4 Risk") +
  scale_color_manual(values = c("Low" = PAL$low, "High" = PAL$high), guide = "none") +
  scale_y_continuous(expand = expansion(mult = c(0.05, 0.18))) +
  labs(x = NULL, y = "ssGSEA Activity (Z-score)") +
  theme_top(8.1) +
  theme(
    legend.position = "bottom", legend.margin = margin(t = -3),
    strip.text = element_text(face = "bold", size = 7.9, margin = margin(t = 4, b = 3)),
    axis.text.x = element_text(size = 7.1),
    plot.margin = margin(12, 6, 4, 6)
  )

# ==============================================================================
# PANEL C: Pan-Checkpoint & Effector Suppression Dumbbell Range
# ==============================================================================
ck <- fread(file.path(res_dir, "GMSig4_ImmuneCheckpoints_3Cohorts.csv"))
ck[, cohort_lbl := factor(cohort,
  levels = c("TCGA_TNBC", "GSE58812", "GSE96058_TNBC_Ext"),
  labels = c("TCGA-TNBC (n=172)", "GSE58812 (n=107)", "GSE96058 (n=242)")
)]
ck_sum <- ck[, .(mean_rho = mean(spearman_rho), min_rho = min(spearman_rho), max_rho = max(spearman_rho)), by = gene][order(mean_rho)]
ck[, gene_f := factor(gene, levels = rev(ck_sum$gene))]
ck_sum[, gene_f := factor(gene, levels = rev(ck_sum$gene))]

p4c <- ggplot() +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_segment(data = ck_sum, aes(x = min_rho, xend = max_rho, y = gene_f, yend = gene_f),
               color = "#94A3B8", linewidth = 1.15, alpha = 0.65) +
  geom_point(data = ck, aes(x = spearman_rho, y = gene_f, color = cohort_lbl, shape = cohort_lbl), size = 2.6) +
  scale_color_manual(values = c(PAL$high, PAL$emer, PAL$navy), name = NULL) +
  scale_shape_manual(values = c(16, 17, 15), name = NULL) +
  scale_x_continuous(limits = c(-0.92, 0.08), breaks = seq(-0.8, 0, 0.2)) +
  scale_y_discrete(labels = parse(text = sprintf("italic('%s')", levels(ck$gene_f)))) +
  labs(x = "Spearman Correlation (\u03C1) w/ GMSig-4", y = NULL) +
  theme_top(8.1) +
  theme(
    axis.text.y = element_text(face = "bold", size = 7.6),
    legend.position = "bottom", legend.margin = margin(t = -3)
  ) +
  guides(color = guide_legend(nrow = 2), shape = guide_legend(nrow = 2))

# ==============================================================================
# PANEL D: TCGA-TNBC Volcano Plot w/ GMSig-4 & Checkpoint Tiered Highlights
# ==============================================================================
deg <- fread(file.path(res_dir, "GMSig4_TCGA_TNBC_HighVsLow_limma_DEGs.csv"))
deg[, neglog10FDR := -log10(pmax(adj.P.Val, 1e-35))]
deg[, Status := ifelse(adj.P.Val < 0.05 & logFC >= 0.75, "Up in High Risk",
                ifelse(adj.P.Val < 0.05 & logFC <= -0.75, "Down in High Risk", "NS"))]
gmsig_genes <- c("IFNG", "IDO1", "CXCL10", "ENO1")
top_down <- deg[gene %in% c("IFNG", "IDO1", "CXCL10", "CXCL9", "CD8A", "PDCD1", "FASLG", "GBP4", "CXCR6", "TBX21", "GZMB")]
top_up   <- deg[adj.P.Val < 0.05 & logFC > 0.75][order(-t)][1:5]
if (!("ENO1" %in% top_up$gene) && "ENO1" %in% deg$gene) top_up <- rbind(top_up, deg[gene == "ENO1"])
lbl_dt <- unique(rbind(top_down, top_up), by = "gene")
lbl_dt[, is_sig4 := gene %in% gmsig_genes]

n_up   <- sum(deg$Status == "Up in High Risk")
n_down <- sum(deg$Status == "Down in High Risk")

p4d <- ggplot(deg, aes(x = logFC, y = neglog10FDR, color = Status)) +
  geom_vline(xintercept = c(-0.75, 0.75), linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "#94A3B8", linewidth = 0.4) +
  geom_point(size = 0.65, alpha = 0.45) +
  geom_point(data = lbl_dt[is_sig4 == FALSE], size = 1.9, shape = 21, fill = PAL$ink, color = "white", stroke = 0.4) +
  geom_point(data = lbl_dt[is_sig4 == TRUE], size = 2.9, shape = 23,
             aes(fill = ifelse(logFC > 0, PAL$high, PAL$low)), color = PAL$ink, stroke = 0.7) +
  scale_fill_identity() +
  ggrepel::geom_text_repel(
    data = lbl_dt, aes(label = gene, fontface = ifelse(is_sig4, "bold.italic", "italic")),
    color = PAL$ink, size = 2.35, family = "Arial",
    max.overlaps = 30, box.padding = 0.28, point.padding = 0.15, min.segment.length = 0.08, segment.size = 0.35
  ) +
  annotate("label", x = -3.6, y = 29, label = sprintf("Down: %d", n_down), hjust = 0,
           size = 2.4, fontface = "bold", family = "Arial", color = PAL$low, fill = "#EEF6FA", linewidth = 0.3) +
  annotate("label", x = 2.3, y = 29, label = sprintf("Up: %d", n_up), hjust = 1,
           size = 2.4, fontface = "bold", family = "Arial", color = PAL$high, fill = "#FDEEEC", linewidth = 0.3) +
  scale_color_manual(values = c("Down in High Risk" = PAL$low, "NS" = "#CBD5E0", "Up in High Risk" = PAL$high), name = NULL) +
  labs(x = expression(bold(log[2]*" Fold Change (High vs Low GMSig-4)")), y = expression(bold(-log[10]*"(BH FDR)"))) +
  theme_top(8.1) +
  theme(legend.position = "bottom", legend.margin = margin(t = -3))

# ==============================================================================
# PANEL E: Multi-Pathway Hallmark GSEA Running Enrichment Curves + Rank Barcode
# ==============================================================================
h_df <- clusterProfiler::read.gmt(file.path(msig_dir, "hallmark.gmt"))
h_sets <- split(h_df$gene, h_df$term)
g_h <- fread(file.path(res_dir, "GMSig4_TCGA_TNBC_GSEA_Hallmarks.csv"))

rank_stats <- sort(setNames(deg$t, deg$gene), decreasing = TRUE)
sel_hm <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION",
  "HALLMARK_GLYCOLYSIS",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"
)
hm_labels <- c(
  "HALLMARK_EPITHELIAL_MESENCHYMAL_TRANSITION" = "EMT (NES=+2.43)",
  "HALLMARK_GLYCOLYSIS"                        = "Glycolysis (NES=+2.03)",
  "HALLMARK_INTERFERON_ALPHA_RESPONSE"         = "IFN-\u03B1 (NES=\u22123.03)",
  "HALLMARK_INTERFERON_GAMMA_RESPONSE"         = "IFN-\u03B3 (NES=\u22123.21)"
)
hm_cols <- c(
  "EMT (NES=+2.43)"        = PAL$high,
  "Glycolysis (NES=+2.03)" = PAL$ochre,
  "IFN-\u03B1 (NES=\u22123.03)" = PAL$low,
  "IFN-\u03B3 (NES=\u22123.21)" = PAL$navy
)

curve_list <- list(); tick_list <- list()
for (i_hm in seq_along(sel_hm)) {
  pw_id <- sel_hm[i_hm]
  pd <- fgsea::plotEnrichmentData(h_sets[[pw_id]], rank_stats)
  c_dt <- as.data.table(pd$curve)
  c_dt[, Pathway := factor(hm_labels[pw_id], levels = unname(hm_labels[sel_hm]))]
  curve_list[[i_hm]] <- c_dt
  t_dt <- as.data.table(pd$ticks)
  t_dt[, Pathway := factor(hm_labels[pw_id], levels = unname(hm_labels[sel_hm]))]
  t_dt[, row_y := length(sel_hm) - i_hm + 1]
  tick_list[[i_hm]] <- t_dt
}
curves_dt <- rbindlist(curve_list)
ticks_dt  <- rbindlist(tick_list)

p4e_top <- ggplot(curves_dt, aes(x = rank, y = ES, color = Pathway)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_line(linewidth = 0.95) +
  scale_color_manual(values = hm_cols, name = NULL) +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, 5000, 10000, 15000, 20000)) +
  labs(x = NULL, y = "Running Enrichment Score") +
  theme_top(8.1) +
  theme(
    axis.text.x = element_blank(), axis.ticks.x = element_blank(),
    legend.position = c(0.75, 0.80),
    legend.background = element_rect(fill = alpha("white", 0.88), color = "#CBD5E0", linewidth = 0.35),
    legend.key.size = unit(9, "pt"),
    plot.margin = margin(6, 6, 0, 6)
  )

p4e_ticks <- ggplot(ticks_dt, aes(x = rank, y = row_y, color = Pathway)) +
  geom_segment(aes(xend = rank, y = row_y - 0.38, yend = row_y + 0.38), linewidth = 0.32, alpha = 0.75) +
  scale_color_manual(values = hm_cols, guide = "none") +
  scale_x_continuous(expand = c(0.01, 0), breaks = c(0, 5000, 10000, 15000, 20000)) +
  scale_y_continuous(breaks = 4:1, labels = c("EMT", "Glycolysis", "IFN-\u03B1", "IFN-\u03B3")) +
  labs(x = "Gene Rank in TCGA-TNBC (High \u2192 Low GMSig-4)", y = NULL) +
  theme_top(7.8) +
  theme(
    axis.line.y = element_blank(), axis.ticks.y = element_blank(),
    axis.text.y = element_text(size = 6.9, face = "bold"),
    plot.margin = margin(1, 6, 4, 6)
  )

p4e <- (p4e_top / p4e_ticks) + plot_layout(heights = c(1.0, 0.35))

# ==============================================================================
# PANEL F: MSigDB KEGG GSEA Diverging Bubble Lollipop
# ==============================================================================
g_k <- fread(file.path(res_dir, "GMSig4_TCGA_TNBC_GSEA_KEGG.csv"))
top_pos <- g_k[p.adjust < 0.05 & NES > 0][order(-abs(NES))][1:min(6, .N)]
top_neg <- g_k[p.adjust < 0.05 & NES < 0][order(-abs(NES))][1:min(8, .N)]
sub_k <- unique(rbind(top_pos, top_neg), by = "ID")
sub_k[, clean := gsub("^KEGG_", "", ID)]
sub_k[, clean := gsub("_", " ", clean)]
sub_k[, clean := tools::toTitleCase(tolower(clean))]
repl <- c(
  "Gpi Anchor" = "GPI-Anchor",
  "Glycosylphosphatidylinositol GPI-Anchor Biosynthesis" = "GPI-Anchor Biosynthesis",
  "Natural Killer Cell Mediated Cytotoxicity" = "NK Cell Mediated Cytotoxicity",
  "Cell Adhesion Molecules Cams" = "Cell Adhesion Molecules (CAMs)",
  "Intestinal Immune Network For" = "Intestinal Immune Network for",
  "Iga " = "IgA ", "Ecm Receptor" = "ECM-Receptor"
)
for (k in names(repl)) sub_k[, clean := gsub(k, repl[[k]], clean)]
sub_k[, Dir := ifelse(NES > 0, "Activated in High Risk", "Suppressed in High Risk")]

p4f <- ggplot(sub_k, aes(x = NES, y = reorder(clean, NES))) +
  geom_vline(xintercept = 0, color = PAL$line, linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = NES, yend = reorder(clean, NES), color = Dir), linewidth = 0.85, alpha = 0.8) +
  geom_point(aes(size = setSize, fill = Dir), shape = 21, color = PAL$ink, stroke = 0.45) +
  geom_text(aes(label = sprintf("%+.2f", NES), hjust = ifelse(NES >= 0, -0.35, 1.35)),
            size = 2.25, fontface = "bold", family = "Arial", color = PAL$ink) +
  scale_fill_manual(values = c("Activated in High Risk" = PAL$high, "Suppressed in High Risk" = PAL$low), name = NULL) +
  scale_color_manual(values = c("Activated in High Risk" = PAL$high, "Suppressed in High Risk" = PAL$low), guide = "none") +
  scale_size_continuous(range = c(2.3, 4.8), breaks = c(50, 100, 150), name = "Genes") +
  scale_x_continuous(limits = c(-3.75, 3.35), breaks = c(-3, -2, -1, 0, 1, 2, 3)) +
  labs(x = "Normalized Enrichment Score (NES)", y = NULL) +
  theme_top(8.1) +
  theme(
    axis.text.y = element_text(size = 7.2),
    legend.position = "bottom", legend.box = "horizontal", legend.margin = margin(t = -3)
  ) +
  guides(fill = guide_legend(nrow = 1, override.aes = list(size = 3.0)), size = guide_legend(nrow = 1))

# ==============================================================================
# Composite Assembly
# ==============================================================================
wrap_panel <- function(plot_obj, tag_letter, title_txt) {
  wrap_elements(full = plot_obj) +
    labs(tag = tag_letter, title = title_txt) +
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 14, color = "#111111"),
      plot.tag.position = c(0.01, 0.986),
      plot.title = element_text(family = "Arial", face = "bold", size = 9.7, color = PAL$ink,
                                hjust = 0, margin = margin(l = 22, t = 2, b = 6)),
      plot.margin = margin(5, 5, 4, 5)
    )
}

w4a <- wrap_panel(p4a, "A", "3-Cohort Consensus TME & Pathway Correlation")
w4b <- wrap_panel(p4b, "B", "High- vs Low-Risk TME & Microbial Axis Depletion")
w4c <- wrap_panel(p4c, "C", "Pan-Checkpoint & Effector Suppression")
w4d <- wrap_panel(p4d, "D", "TCGA-TNBC High vs Low GMSig-4 Transcriptome")
w4e <- wrap_panel(p4e, "E", "Hallmark GSEA Enrichment Curves (TCGA-TNBC)")
w4f <- wrap_panel(p4f, "F", "MSigDB KEGG GSEA Pathways (TCGA-TNBC)")

fig4 <- (w4a | w4b) / (w4c | w4d) / (w4e | w4f) +
  plot_layout(heights = c(1.12, 1.05, 1.05), widths = c(1.05, 1.05))

agg_png(file.path(fig_dir, "Figure4.png"), width = 15.8, height = 15.6, units = "in", res = 300, background = "white")
print(fig4)
dev.off()

cairo_pdf(file.path(fig_dir, "Figure4.pdf"), width = 15.8, height = 15.6)
print(fig4)
dev.off()

cat("SAVED Figure4.png and Figure4.pdf successfully!\n")
