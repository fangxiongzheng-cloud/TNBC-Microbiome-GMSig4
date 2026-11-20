project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
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

res_dir <- file.path(project_root, "02_results")
fig_dir <- file.path(project_root, "03_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

PAL <- list(
  overall = "#64748B",
  erneg   = "#2B7A9E",
  tnbc    = "#C84B31",
  pure    = "#8B1E1E",
  emer    = "#1F8A70",
  navy    = "#1E3A5F",
  div_lo  = "#1D5B9B",
  div_mi  = "#F8FAFC",
  div_hi  = "#C0392B",
  ink     = "#1A202C",
  muted   = "#4A5568",
  line    = "#2D3748"
)

theme_top <- function(base_size = 8.3) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", color = PAL$ink),
      axis.line = element_line(color = PAL$line, linewidth = 0.45),
      axis.ticks = element_line(color = PAL$line, linewidth = 0.4),
      axis.ticks.length = unit(2.4, "pt"),
      axis.title = element_text(face = "bold", size = base_size, color = PAL$ink),
      axis.text = element_text(size = base_size - 0.8, color = "#2D3748"),
      plot.title = element_text(face = "bold", size = base_size + 1.3, hjust = 0, margin = margin(b = 5)),
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
# Load 6-GWAS Combined MR Tables
# ==============================================================================
gwas_levels <- c("ieu-a-1126", "GCST010098", "ieu-a-1128", "ieu-a-1135", "GCST010100", "GCST90454344")
gwas_labels <- c(
  "ieu-a-1126"   = "Overall BC\n(BCAC 1126)",
  "GCST010098"   = "Overall BC\n(Zhang)",
  "ieu-a-1128"   = "ER\u2212 BC\n(Meta 1128)",
  "ieu-a-1135"   = "ER\u2212 BC\n(Onco 1135)",
  "GCST010100"   = "TNBC\n(BCAC)",
  "GCST90454344" = "Pure TNBC\n(90454344)"
)
subtype_grp <- c(
  "ieu-a-1126"   = "Overall BC",
  "GCST010098"   = "Overall BC",
  "ieu-a-1128"   = "ER\u2212 BC",
  "ieu-a-1135"   = "ER\u2212 BC",
  "GCST010100"   = "TNBC",
  "GCST90454344" = "TNBC"
)

mib5  <- fread(file.path(res_dir, "MR_MiBioGen_5GWAS_CrossSubtype_Long.csv"))
mib98 <- fread(file.path(res_dir, "MR_MiBioGen_GCST010098_results.csv"))[
  method %in% c("Inverse variance weighted", "Wald ratio"),
  .(outcome_id = id.outcome, exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)
]
mib6 <- rbind(mib5, mib98)

eq5  <- fread(file.path(res_dir, "MR_eQTLGen_5GWAS_CrossSubtype_Long.csv"))
eq98 <- fread(file.path(res_dir, "MR_eQTLGen_GCST010098_results.csv"))[
  method %in% c("Inverse variance weighted", "Wald ratio"),
  .(outcome_id = id.outcome, exposure, nsnp, OR, OR_lci, OR_uci, pval, fdr)
]
eq6 <- rbind(eq5, eq98)

clean_taxon_expr <- function(x) {
  y <- gsub("\\.id\\.[0-9]+$", "", x)
  y <- gsub("^genus\\.\\.?", "g. ", y)
  y <- gsub("^family\\.", "f. ", y)
  y <- gsub("^order\\.", "o. ", y)
  y <- gsub("Eubacteriumxylanophilumgroup", "Eubact. xylanophilum", y)
  y <- gsub("LachnospiraceaeNK4A136group", "Lachno. NK4A136", y)
  y <- gsub("FamilyXIIIUCG001", "Family XIII UCG001", y)
  y
}

# ==============================================================================
# PANEL A: Annotated Cross-Subtype Microbiota Causal Matrix
# ==============================================================================
sel_taxa <- c(
  "genus.Dorea.id.1997",
  "genus.FamilyXIIIUCG001.id.11294",
  "genus.Desulfovibrio.id.3173",
  "genus.Parabacteroides.id.954",
  "genus.Marvinbryantia.id.2005",
  "genus.Oxalobacter.id.2978",
  "genus.LachnospiraceaeNK4A136group.id.11319",
  "genus..Eubacteriumxylanophilumgroup.id.14375",
  "order.Lactobacillales.id.1800",
  "genus.Adlercreutzia.id.812",
  "genus.Sellimonas.id.14369",
  "family.Ruminococcaceae.id.2050",
  "genus.Anaerostipes.id.1991"
)
m2a <- mib6[exposure %in% sel_taxa]
m2a[, taxon_plain := clean_taxon_expr(exposure)]
taxon_ord <- rev(clean_taxon_expr(sel_taxa))
m2a[, taxon_f := factor(taxon_plain, levels = taxon_ord)]
m2a[, gwas_f := factor(outcome_id, levels = gwas_levels, labels = gwas_labels[gwas_levels])]
m2a[, log2OR := log2(OR)]
m2a[, cell_txt := sprintf("%.2f%s", OR, ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", "")))]
m2a[, txt_col := ifelse(abs(log2OR) > 0.25, "white", PAL$ink)]

top_strip_a <- data.table(
  gwas_f = factor(gwas_levels, levels = gwas_levels, labels = gwas_labels[gwas_levels]),
  Subtype = factor(subtype_grp[gwas_levels], levels = c("Overall BC", "ER\u2212 BC", "TNBC"))
)
p2a_top <- ggplot(top_strip_a, aes(x = gwas_f, y = 1, fill = Subtype)) +
  geom_tile(color = "white", linewidth = 0.8, height = 0.8) +
  geom_text(aes(label = Subtype), color = "white", fontface = "bold", size = 2.3, family = "Arial") +
  scale_fill_manual(values = c("Overall BC" = PAL$overall, "ER\u2212 BC" = PAL$erneg, "TNBC" = PAL$tnbc), guide = "none") +
  theme_void() +
  theme(plot.margin = margin(8, 6, 0, 6))

italic_taxa_a <- parse(text = sapply(taxon_ord, function(s) {
  sub("^([^.]+\\.)\\s*(.*)$", "'\\1'~italic('\\2')", s)
}))

p2a_mat <- ggplot(m2a, aes(x = gwas_f, y = taxon_f, fill = log2OR)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_tile(data = m2a[pval < 0.05], fill = NA, color = PAL$ink, linewidth = 0.65) +
  geom_text(aes(label = cell_txt, color = txt_col, fontface = ifelse(pval < 0.05, "bold", "plain")),
            size = 2.4, family = "Arial") +
  scale_color_identity() +
  scale_fill_gradient2(low = PAL$div_lo, mid = PAL$div_mi, high = PAL$div_hi, midpoint = 0,
                       limits = c(-0.45, 0.35), name = expression(bold(log[2]*"(OR)"))) +
  scale_y_discrete(labels = italic_taxa_a) +
  labs(x = NULL, y = NULL) +
  theme_top(8.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(size = 6.9, face = "bold"),
    axis.text.y = element_text(size = 7.4),
    legend.position = "bottom", legend.margin = margin(t = -3),
    plot.margin = margin(1, 6, 4, 6)
  ) +
  guides(fill = guide_colorbar(barwidth = unit(100, "pt"), barheight = unit(5.5, "pt"),
                               title.position = "left", title.vjust = 0.85))

p2a <- (p2a_top / p2a_mat) + plot_layout(heights = c(0.08, 1.0))

# ==============================================================================
# PANEL B: Subtype Effect Gradient Forest Plot (Zebra Shaded)
# ==============================================================================
grad_taxa <- c(
  "genus.Dorea.id.1997",
  "genus.Sellimonas.id.14369",
  "genus.FamilyXIIIUCG001.id.11294",
  "genus.Desulfovibrio.id.3173",
  "genus.Marvinbryantia.id.2005",
  "genus.Anaerostipes.id.1991"
)
sub_gwas <- c("ieu-a-1126", "ieu-a-1128", "GCST010100", "GCST90454344")
sub_glbl <- c(
  "ieu-a-1126"   = "Overall BC (1126)",
  "ieu-a-1128"   = "ER\u2212 BC (1128)",
  "GCST010100"   = "TNBC (010100)",
  "GCST90454344" = "Pure TNBC (90454344)"
)
f2b <- mib6[exposure %in% grad_taxa & outcome_id %in% sub_gwas]
ord_b <- rev(clean_taxon_expr(grad_taxa))
f2b[, taxon_f := factor(clean_taxon_expr(exposure), levels = ord_b)]
f2b[, Subtype := factor(outcome_id, levels = rev(sub_gwas), labels = rev(sub_glbl[sub_gwas]))]

zebra_b <- data.table(y = seq_along(ord_b))[y %% 2 == 1]
italic_taxa_b <- parse(text = sapply(ord_b, function(s) {
  sub("^([^.]+\\.)\\s*(.*)$", "'\\1'~italic('\\2')", s)
}))

p2b <- ggplot(f2b, aes(x = OR, y = taxon_f, color = Subtype)) +
  geom_rect(data = zebra_b, aes(ymin = y - 0.48, ymax = y + 0.48, xmin = 0.56, xmax = 1.50),
            inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_errorbar(aes(xmin = OR_lci, xmax = OR_uci), orientation = "y", width = 0.22, linewidth = 0.7,
                position = position_dodge(width = 0.68)) +
  geom_point(size = 2.1, position = position_dodge(width = 0.68)) +
  scale_color_manual(
    values = c("Overall BC (1126)" = PAL$overall, "ER\u2212 BC (1128)" = PAL$erneg,
               "TNBC (010100)" = PAL$tnbc, "Pure TNBC (90454344)" = PAL$pure),
    breaks = unname(sub_glbl[sub_gwas]), name = NULL
  ) +
  scale_x_log10(breaks = c(0.65, 0.8, 1.0, 1.2, 1.45), limits = c(0.56, 1.50)) +
  scale_y_discrete(labels = italic_taxa_b) +
  labs(x = "Mendelian Randomization OR (95% CI)", y = NULL) +
  theme_top(8.2) +
  theme(
    legend.position = "bottom", legend.key.size = unit(8, "pt"),
    legend.margin = margin(t = -3), axis.text.y = element_text(size = 7.6)
  ) +
  guides(color = guide_legend(nrow = 2))

# ==============================================================================
# PANEL C: Host gutMGene2 Causal Targets (eQTLGen -> 6 BC GWAS)
# ==============================================================================
sel_genes <- c("PEMT", "CCR6", "MUC1", "TCF7L2", "FFAR2", "NLRP6", "ICAM1", "FOXP3", "MAGI3", "ST6GAL1", "CXCR6")
gene_cls <- c(
  "PEMT" = "Lipid / Metabolic", "FFAR2" = "Lipid / Metabolic", "TCF7L2" = "Lipid / Metabolic", "ST6GAL1" = "Lipid / Metabolic",
  "CCR6" = "Immune / Barrier", "NLRP6" = "Immune / Barrier", "FOXP3" = "Immune / Barrier", "CXCR6" = "Immune / Barrier",
  "MUC1" = "Epithelial / Junction", "ICAM1" = "Epithelial / Junction", "MAGI3" = "Epithelial / Junction"
)
gcls_cols <- c("Lipid / Metabolic" = PAL$emer, "Immune / Barrier" = PAL$erneg, "Epithelial / Junction" = "#8E6C8A")

e2c <- eq6[exposure %in% sel_genes]
e2c[, gene_f := factor(exposure, levels = rev(sel_genes))]
e2c[, gwas_f := factor(outcome_id, levels = gwas_levels, labels = gwas_labels[gwas_levels])]
e2c[, log2OR := pmax(pmin(log2(OR), 0.75), -0.60)]
e2c[, cell_txt := sprintf("%.2f%s", OR,
                          ifelse(!is.na(fdr) & fdr < 0.05, "***",
                          ifelse(pval < 0.01, "**", ifelse(pval < 0.05, "*", ""))))]
e2c[, txt_col := ifelse(abs(log2OR) > 0.36, "white", PAL$ink)]

left_strip_c <- data.table(
  gene_f = factor(rev(sel_genes), levels = rev(sel_genes)),
  Module = factor(gene_cls[rev(sel_genes)], levels = names(gcls_cols))
)

p2c_top <- ggplot(top_strip_a, aes(x = gwas_f, y = 1, fill = Subtype)) +
  geom_tile(color = "white", linewidth = 0.8, height = 0.8) +
  geom_text(aes(label = Subtype), color = "white", fontface = "bold", size = 2.3, family = "Arial") +
  scale_fill_manual(values = c("Overall BC" = PAL$overall, "ER\u2212 BC" = PAL$erneg, "TNBC" = PAL$tnbc), guide = "none") +
  theme_void() +
  theme(plot.margin = margin(8, 6, 0, 1))

p2c_left <- ggplot(left_strip_c, aes(x = 1, y = gene_f, fill = Module)) +
  geom_tile(color = "white", linewidth = 0.6, width = 0.85) +
  scale_fill_manual(values = gcls_cols, name = NULL) +
  scale_y_discrete(labels = parse(text = sprintf("italic('%s')", rev(sel_genes)))) +
  labs(x = NULL, y = NULL) +
  theme_top(8.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7.5, face = "bold"),
    legend.position = "bottom", legend.key.size = unit(7, "pt"),
    plot.margin = margin(1, 0, 4, 4)
  ) +
  guides(fill = guide_legend(nrow = 2))

p2c_mat <- ggplot(e2c, aes(x = gwas_f, y = gene_f, fill = log2OR)) +
  geom_tile(color = "white", linewidth = 0.7) +
  geom_tile(data = e2c[pval < 0.05], fill = NA, color = PAL$ink, linewidth = 0.65) +
  geom_text(aes(label = cell_txt, color = txt_col, fontface = ifelse(pval < 0.05, "bold", "plain")),
            size = 2.35, family = "Arial") +
  scale_color_identity() +
  scale_fill_gradient2(low = PAL$div_lo, mid = PAL$div_mi, high = PAL$div_hi, midpoint = 0,
                       limits = c(-0.62, 0.78), name = expression(bold(log[2]*"(OR)"))) +
  labs(x = NULL, y = NULL) +
  theme_top(8.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.y = element_blank(), axis.text.x = element_text(size = 6.9, face = "bold"),
    legend.position = "bottom", plot.margin = margin(1, 6, 4, 1)
  ) +
  guides(fill = guide_colorbar(barwidth = unit(72, "pt"), barheight = unit(5.5, "pt"),
                               title.position = "left", title.vjust = 0.85))

p2c <- wrap_plots(plot_spacer(), p2c_top, p2c_left, p2c_mat,
                  ncol = 2, nrow = 2, heights = c(0.08, 1.0), widths = c(0.11, 1.0), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(6, "pt"))

# ==============================================================================
# PANEL D: GTEx V8 Breast Colocalization Dumbbell (Marginal PP.H4 vs cond_H4)
# ==============================================================================
coloc <- fread(file.path(res_dir, "Coloc_GTExBreast_TNBC_results.csv"))
gt_mr <- fread(file.path(res_dir, "MR_GTExBreast_GCST010100_TNBC_results.csv"))[method %in% c("Inverse variance weighted", "Wald ratio")]
d2 <- merge(coloc, gt_mr[, .(gene = exposure, GTEx_OR = OR, GTEx_p = pval)], by = "gene", all.x = TRUE)
d2_sub <- d2[gene %in% c("TCF7L2", "CLDN7", "AKT1", "ENO1", "PEMT", "FFAR2", "VEGFA", "SFN", "MUC1", "DAPK3", "ICAM1", "LYZ")]
d2_sub[, top_col := cond_H4 >= 0.35]
d2_sub <- d2_sub[order(cond_H4)]
d2_sub[, gene_f := factor(gene, levels = gene)]
d2_sub[, lbl_r := sprintf("cond %.0f%% | H4 %.1f%%", cond_H4 * 100, PP.H4 * 100)]

p2d <- ggplot(d2_sub, aes(y = gene_f)) +
  geom_vline(xintercept = 35, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_segment(aes(x = PP.H4 * 100, xend = cond_H4 * 100, yend = gene_f, color = top_col), linewidth = 1.1) +
  geom_point(aes(x = PP.H4 * 100, shape = "Marginal PP.H4"), size = 2.4, fill = "white", color = PAL$navy, stroke = 0.8) +
  geom_point(aes(x = cond_H4 * 100, fill = top_col, shape = "Conditional H4"), size = 3.1, color = PAL$ink, stroke = 0.4) +
  geom_text(aes(x = cond_H4 * 100 + 2.5, label = lbl_r), hjust = 0, size = 2.2, fontface = "bold", family = "Arial", color = PAL$ink) +
  scale_color_manual(values = c("FALSE" = "#94A3B8", "TRUE" = PAL$emer), guide = "none") +
  scale_fill_manual(values = c("FALSE" = "#94A3B8", "TRUE" = PAL$emer),
                    labels = c("FALSE" = "cond_H4 < 35%", "TRUE" = "cond_H4 \u2265 35%"), name = NULL) +
  scale_shape_manual(values = c("Marginal PP.H4" = 21, "Conditional H4" = 21), name = NULL) +
  scale_x_continuous(limits = c(-2, 114), breaks = c(0, 25, 50, 75, 100)) +
  scale_y_discrete(labels = parse(text = sprintf("italic('%s')", levels(d2_sub$gene_f)))) +
  labs(x = "Colocalization Posterior Probability (%)", y = NULL) +
  theme_top(8.2) +
  theme(
    legend.position = "bottom", legend.margin = margin(t = -3),
    axis.text.y = element_text(face = "bold", size = 7.6)
  ) +
  guides(
    shape = guide_legend(override.aes = list(fill = c("white", PAL$emer), size = 2.6)),
    fill = guide_legend(override.aes = list(shape = 22, size = 3.1))
  )

# ==============================================================================
# PANEL E: GTEx V8 Breast cis-eQTL MR -> TNBC Forest Table Plot
# ==============================================================================
gt_sub <- gt_mr[order(pval)][1:min(14, .N)]
gt_sub <- gt_sub[order(-pval)]
gt_sub[, gene_f := factor(exposure, levels = exposure)]
gt_sub[, sig_col := factor(
  ifelse(pval < 0.01, "P < 0.01 (FDR=0.083)", ifelse(pval < 0.10, "P < 0.10", "NS")),
  levels = c("P < 0.01 (FDR=0.083)", "P < 0.10", "NS")
)]
gt_sub[, tbl_or := sprintf("%.2f (%.2f\u2013%.2f)", OR, OR_lci, OR_uci)]
gt_sub[, tbl_p  := ifelse(pval < 0.01, sprintf("%.4f", pval), sprintf("%.3f", pval))]
zebra_e <- data.table(y = seq_len(nrow(gt_sub)))[y %% 2 == 1]

p2e_left <- ggplot(gt_sub, aes(y = gene_f)) +
  geom_text(aes(x = 1, label = nsnp), size = 2.4, family = "Arial", color = PAL$muted) +
  scale_y_discrete(labels = parse(text = sprintf("italic('%s')", levels(gt_sub$gene_f)))) +
  scale_x_continuous(limits = c(0.6, 1.4), breaks = 1, labels = "IVs") +
  labs(x = NULL, y = NULL) +
  theme_top(8.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 7.2, color = PAL$muted),
    axis.text.y = element_text(face = "bold", size = 7.5),
    plot.margin = margin(6, 0, 4, 4)
  )

p2e_mid <- ggplot(gt_sub, aes(x = OR, y = gene_f, color = sig_col)) +
  geom_rect(data = zebra_e, aes(ymin = y - 0.48, ymax = y + 0.48, xmin = 0.55, xmax = 4.6),
            inherit.aes = FALSE, fill = "#F8FAFC", color = NA) +
  geom_vline(xintercept = 1, linetype = "dashed", color = "#718096", linewidth = 0.45) +
  geom_errorbar(aes(xmin = OR_lci, xmax = OR_uci), orientation = "y", width = 0.22, linewidth = 0.7) +
  geom_point(size = 2.4) +
  scale_x_log10(breaks = c(0.6, 0.8, 1.0, 1.5, 2.5, 4.0), limits = c(0.55, 4.6)) +
  scale_color_manual(values = c("P < 0.01 (FDR=0.083)" = PAL$tnbc, "P < 0.10" = PAL$emer, "NS" = PAL$erneg), name = NULL) +
  labs(x = "OR (95% CI)", y = NULL) +
  theme_top(8.2) +
  theme(
    axis.text.y = element_blank(), axis.ticks.y = element_blank(), axis.line.y = element_blank(),
    legend.position = "bottom", legend.margin = margin(t = -3),
    plot.margin = margin(6, 2, 4, 0)
  )

p2e_tbl <- ggplot(gt_sub, aes(y = gene_f)) +
  geom_text(aes(x = 0.0, label = tbl_or, fontface = ifelse(pval < 0.10, "bold", "plain")),
            hjust = 0, size = 2.35, family = "Arial", color = PAL$ink) +
  geom_text(aes(x = 1.15, label = tbl_p, fontface = ifelse(pval < 0.10, "bold", "plain"),
                color = ifelse(pval < 0.01, PAL$tnbc, PAL$ink)),
            hjust = 1, size = 2.35, family = "Arial") +
  scale_color_identity() +
  scale_x_continuous(limits = c(-0.02, 1.2), breaks = c(0.4, 1.05), labels = c("OR (95% CI)", "P")) +
  labs(x = NULL, y = NULL) +
  theme_top(8.2) +
  theme(
    axis.line = element_blank(), axis.ticks = element_blank(), axis.text.y = element_blank(),
    axis.text.x = element_text(face = "bold", size = 7.2, color = PAL$muted),
    plot.margin = margin(6, 4, 4, 0)
  )

p2e <- (p2e_left | p2e_mid | p2e_tbl) + plot_layout(widths = c(0.24, 1.22, 0.94))

# ==============================================================================
# PANEL F: Two-Step Mediation MR (3-Node Causal Path DAG + Effect Bars)
# ==============================================================================
med <- fread(file.path(res_dir, "MR_Mediation_2Step_Causal_Pathways.csv"))

build_med_panel <- function(row_i) {
  r <- med[row_i]
  tax_lbl <- clean_taxon_expr(r$taxon)
  oc_lbl  <- ifelse(r$outcome == "GCST010100", "TNBC (010100)", "Overall BC (1126)")

  p_dag <- ggplot() +
    xlim(0, 10) + ylim(0, 3.6) +
    annotate("rect", xmin = 0.1, xmax = 3.1, ymin = 1.1, ymax = 2.4, fill = "#EEF6FA", color = PAL$erneg, linewidth = 0.65) +
    annotate("text", x = 1.6, y = 1.75, label = tax_lbl, fontface = "bold.italic", size = 2.35, family = "Arial", color = PAL$ink) +
    annotate("rect", xmin = 3.8, xmax = 6.2, ymin = 1.1, ymax = 2.4, fill = "#ECF7F4", color = PAL$emer, linewidth = 0.65) +
    annotate("text", x = 5.0, y = 1.75, label = paste0(r$mediator_gene, " eQTL"), fontface = "bold", size = 2.35, family = "Arial", color = PAL$ink) +
    annotate("rect", xmin = 6.9, xmax = 9.9, ymin = 1.1, ymax = 2.4, fill = "#FDEEEC", color = PAL$tnbc, linewidth = 0.65) +
    annotate("text", x = 8.4, y = 1.75, label = oc_lbl, fontface = "bold", size = 2.35, family = "Arial", color = PAL$ink) +
    annotate("segment", x = 3.15, xend = 3.75, y = 1.75, yend = 1.75,
             arrow = arrow(length = unit(4, "pt"), type = "closed"), color = PAL$ink, linewidth = 0.55) +
    annotate("segment", x = 6.25, xend = 6.85, y = 1.75, yend = 1.75,
             arrow = arrow(length = unit(4, "pt"), type = "closed"), color = PAL$ink, linewidth = 0.55) +
    annotate("text", x = 3.45, y = 2.82, label = sprintf("b1=%+.2f (P=%.3f)", r$b_step1, r$p_step1),
             size = 2.1, fontface = "bold", family = "Arial", color = PAL$navy) +
    annotate("text", x = 6.55, y = 2.82, label = sprintf("b2=%+.2f (P=%.3f)", r$b_step2, r$p_step2),
             size = 2.1, fontface = "bold", family = "Arial", color = PAL$navy) +
    annotate("segment", x = 1.6, xend = 8.4, y = 0.62, yend = 0.62, linetype = "dashed", color = PAL$muted, linewidth = 0.45) +
    annotate("text", x = 5.0, y = 0.28,
             label = sprintf("Indirect b = %+.3f (Sobel P = %.3f) | Total b = %+.3f (P = %.3f)",
                             r$b_indirect, r$Sobel_p, r$b_total, r$p_total),
             size = 2.1, fontface = "bold", family = "Arial", color = PAL$muted) +
    theme_void() +
    theme(plot.margin = margin(8, 2, 0, 2))

  dt <- data.table(
    comp = factor(c("Total", "Indirect", "Step 2", "Step 1"), levels = c("Total", "Indirect", "Step 2", "Step 1")),
    beta = c(r$b_total, r$b_indirect, r$b_step2, r$b_step1)
  )
  dt[, Dir := ifelse(beta >= 0, "Risk (+)", "Protective (\u2212)")]
  p_bar <- ggplot(dt, aes(x = beta, y = comp, fill = Dir)) +
    geom_vline(xintercept = 0, color = PAL$line, linewidth = 0.4) +
    geom_col(width = 0.58, color = "white") +
    geom_text(aes(label = sprintf("%+.3f", beta), hjust = ifelse(beta >= 0, -0.15, 1.15)),
              size = 2.2, fontface = "bold", family = "Arial", color = PAL$ink) +
    scale_fill_manual(values = c("Protective (\u2212)" = PAL$div_lo, "Risk (+)" = PAL$div_hi), guide = "none") +
    scale_x_continuous(limits = c(-0.85, 0.65)) +
    labs(x = if (row_i == 1) "Causal Effect (\u03B2)" else NULL, y = NULL) +
    theme_top(7.8) +
    theme(plot.margin = margin(0, 6, 4, 6))

  wrap_plots(p_dag, p_bar, ncol = 1, heights = c(0.85, 1.0))
}

p2f <- wrap_plots(build_med_panel(2), build_med_panel(1), ncol = 1)

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

w2a <- wrap_panel(p2a, "A", "Cross-Subtype Microbiota Causal MR (MiBioGen \u2192 6 BC GWAS)")
w2b <- wrap_panel(p2b, "B", "Subtype Effect Gradient (Overall \u2192 ER\u2212 \u2192 TNBC)")
w2c <- wrap_panel(p2c, "C", "Host gutMGene2 Causal Targets (eQTLGen cis-eQTL)")
w2d <- wrap_panel(p2d, "D", "GTEx V8 Breast Tissue Colocalization (N = 396)")
w2e <- wrap_panel(p2e, "E", "GTEx Breast cis-eQTL MR \u2192 TNBC (GCST010100)")
w2f <- wrap_panel(p2f, "F", "Suggestive Two-Step Mediation MR Circuits")

fig2 <- (w2a | w2b) / (w2c | w2d) / (w2e | w2f) +
  plot_layout(heights = c(1.06, 1.10, 1.08), widths = c(1.06, 1.0))

agg_png(file.path(fig_dir, "Figure2.png"), width = 15.6, height = 16.0, units = "in", res = 300, background = "white")
print(fig2)
dev.off()

cairo_pdf(file.path(fig_dir, "Figure2.pdf"), width = 15.6, height = 16.0)
print(fig2)
dev.off()

cat("SAVED Figure2.png and Figure2.pdf successfully!\n")
