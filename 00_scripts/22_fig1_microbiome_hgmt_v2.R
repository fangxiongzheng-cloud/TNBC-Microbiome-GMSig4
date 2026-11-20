project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
  library(ragg)
  library(systemfonts)
  library(scales)
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
hgmt_dir <- file.path(project_root, "01_data", "HGMT")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

PAL <- list(
  bc     = "#C84B31",
  benign = "#E8A87C",
  hc     = "#2B7A9E",
  other  = "#CBD5E0",
  navy   = "#1E3A5F",
  teal   = "#1F8A70",
  div_lo = "#1D5B9B",
  div_mi = "#F8FAFC",
  div_hi = "#C0392B",
  ink    = "#1A202C",
  muted  = "#4A5568",
  line   = "#2D3748"
)

theme_top <- function(base_size = 8.4) {
  theme_classic(base_size = base_size, base_family = "Arial") +
    theme(
      text = element_text(family = "Arial", color = PAL$ink),
      axis.line = element_line(color = PAL$line, linewidth = 0.45),
      axis.ticks = element_line(color = PAL$line, linewidth = 0.4),
      axis.ticks.length = unit(2.5, "pt"),
      axis.title = element_text(face = "bold", size = base_size, color = PAL$ink),
      axis.text = element_text(size = base_size - 0.8, color = "#2D3748"),
      plot.title = element_text(face = "bold", size = base_size + 1.4, hjust = 0, margin = margin(b = 5)),
      plot.title.position = "plot",
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
# PANEL A: Dual-Ring Sunburst Donut + Aligned Cohort Metadata & Sample Tracks
# ==============================================================================
meta8 <- fread(file.path(res_dir, "HGMT_all8_cohorts_metadata_summary.csv"))
meta8[is.na(country) | country == "", country := "NR"]
meta8 <- meta8[order(n_runs_bacteria)]
ord_levels <- meta8$project_id
meta8[, project_f := factor(project_id, levels = ord_levels)]
meta8[, Assay := factor(assay_type, levels = c("WGS", "16S"))]

# A0: Nested Sunburst Donut (Inner = Phenotype, Outer = Assay)
tot_runs <- sum(meta8$n_runs_bacteria)
tot_bc   <- sum(meta8$n_BC)
tot_hc   <- sum(meta8$n_Health)
tot_ben  <- sum(meta8$n_BreastDiseases)
tot_oth  <- tot_runs - tot_bc - tot_hc - tot_ben

inner_df <- data.table(
  ring = "Phenotype",
  cat  = c("Healthy", "Benign", "BC", "Other"),
  val  = c(tot_hc, tot_ben, tot_bc, tot_oth),
  col  = c(PAL$hc, PAL$benign, PAL$bc, PAL$other)
)
inner_df[, ymax := cumsum(val)][, ymin := shift(ymax, fill = 0)]
inner_df[, pct := val / tot_runs * 100]
inner_df[, ymid := (ymin + ymax) / 2]

tot_16s <- sum(meta8[assay_type == "16S", n_runs_bacteria])
tot_wgs <- sum(meta8[assay_type == "WGS", n_runs_bacteria])
outer_df <- data.table(
  ring = "Assay",
  cat  = c("16S rRNA", "WGS"),
  val  = c(tot_16s, tot_wgs),
  col  = c("#4A7C9B", PAL$navy)
)
outer_df[, ymax := cumsum(val)][, ymin := shift(ymax, fill = 0)]
outer_df[, pct := val / tot_runs * 100]
outer_df[, ymid := (ymin + ymax) / 2]

pA_donut <- ggplot() +
  geom_rect(data = inner_df, aes(xmin = 2.05, xmax = 3.12, ymin = ymin, ymax = ymax),
            fill = inner_df$col, color = "white", linewidth = 0.9) +
  geom_rect(data = outer_df, aes(xmin = 3.22, xmax = 4.22, ymin = ymin, ymax = ymax),
            fill = outer_df$col, color = "white", linewidth = 0.9) +
  geom_text(data = inner_df[pct > 6], aes(x = 2.58, y = ymid, label = sprintf("%.0f%%", pct)),
            color = "white", fontface = "bold", size = 2.25, family = "Arial") +
  geom_text(data = outer_df, aes(x = 3.72, y = ymid, label = sprintf("%.0f%%", pct)),
            color = "white", fontface = "bold", size = 2.25, family = "Arial") +
  annotate("text", x = 0.0, y = 0, label = sprintf("N = %s\nRuns", comma(tot_runs)),
           fontface = "bold", size = 2.7, family = "Arial", color = PAL$ink, lineheight = 1.05) +
  annotate("text", x = 4.9, y = tot_runs * 0.5, label = "", size = 1) +
  coord_polar(theta = "y", start = 0.4) +
  xlim(0.0, 4.35) +
  labs(caption = "Outer: 16S (85%) vs WGS (15%)") +
  theme_void(base_family = "Arial") +
  theme(
    plot.caption = element_text(size = 6.5, face = "bold", color = PAL$navy, hjust = 0.5, margin = margin(t = -6)),
    plot.margin = margin(2, 4, 2, 2)
  )

# A1: Metadata pills
pA_meta <- ggplot(meta8, aes(y = project_f)) +
  geom_tile(aes(x = 1, fill = Assay), width = 0.84, height = 0.66, color = "white", linewidth = 0.5) +
  geom_text(aes(x = 1, label = assay_type), color = "white", fontface = "bold", size = 2.25, family = "Arial") +
  geom_tile(aes(x = 2.02), fill = "#EDF2F7", width = 0.98, height = 0.66, color = "white", linewidth = 0.5) +
  geom_text(aes(x = 2.02, label = country), color = PAL$ink, size = 2.25, family = "Arial") +
  geom_point(data = meta8[has_KO == TRUE], aes(x = 2.82), shape = 18, size = 2.6, color = PAL$teal) +
  geom_point(data = meta8[has_KO == FALSE], aes(x = 2.82), shape = 1, size = 1.8, color = "#CBD5E0") +
  scale_fill_manual(values = c("WGS" = PAL$navy, "16S" = "#4A7C9B"), guide = "none") +
  scale_x_continuous(breaks = c(1, 2.02, 2.82), labels = c("Assay", "Origin", "KO"), expand = c(0.04, 0.04)) +
  labs(x = NULL, y = NULL) +
  theme_top(8.0) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_text(face = "bold", size = 7.4, color = PAL$ink),
    axis.text.x = element_text(face = "bold", size = 7.0, color = PAL$muted),
    plot.margin = margin(4, 6, 4, 2)
  )

# A2: 100% Proportional Stacked Bar
m8_long <- melt(meta8[, .(project_f, n_runs_bacteria,
                          `Breast Cancer` = n_BC,
                          `Benign` = n_BreastDiseases,
                          `Healthy` = n_Health)],
                id.vars = c("project_f", "n_runs_bacteria"), variable.name = "Group", value.name = "Count")
m8_long <- m8_long[Count > 0]
m8_long[, Prop := Count / n_runs_bacteria * 100]
other_494 <- meta8[project_id == "PRJNA494824", n_runs_bacteria - (n_BC + n_Health + n_BreastDiseases)]
if (other_494 > 0) {
  m8_long <- rbind(m8_long, data.table(
    project_f = factor("PRJNA494824", levels = ord_levels),
    n_runs_bacteria = 31, Group = "Other", Count = other_494, Prop = other_494 / 31 * 100
  ))
}
m8_long[, Group := factor(Group, levels = c("Other", "Breast Cancer", "Benign", "Healthy"))]

pA_prop <- ggplot(m8_long, aes(x = Prop, y = project_f, fill = Group)) +
  geom_col(width = 0.66, color = "white", linewidth = 0.4) +
  scale_fill_manual(
    values = c("Healthy" = PAL$hc, "Benign" = PAL$benign, "Breast Cancer" = PAL$bc, "Other" = PAL$other),
    breaks = c("Healthy", "Benign", "Breast Cancer", "Other"),
    name = NULL
  ) +
  scale_x_continuous(breaks = c(0, 50, 100), labels = c("0%", "50%", "100%"), expand = c(0, 0)) +
  labs(x = "Composition", y = NULL) +
  theme_top(8.0) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.key.size = unit(6.5, "pt"),
    legend.margin = margin(t = -3),
    plot.margin = margin(4, 2, 4, 1)
  ) +
  guides(fill = guide_legend(nrow = 2))

# A3: Log10 Sample Lollipop
pA_n <- ggplot(meta8, aes(x = n_runs_bacteria, y = project_f)) +
  geom_segment(aes(x = 15, xend = n_runs_bacteria, yend = project_f), color = "#94A3B8", linewidth = 0.85) +
  geom_point(aes(color = Assay), size = 2.8) +
  geom_text(aes(label = comma(n_runs_bacteria)), hjust = -0.32, size = 2.4, fontface = "bold", family = "Arial") +
  scale_color_manual(values = c("WGS" = PAL$navy, "16S" = "#4A7C9B"), guide = "none") +
  scale_x_log10(limits = c(15, 4200), breaks = c(30, 100, 300, 1000)) +
  labs(x = "Runs (log10)", y = NULL) +
  theme_top(8.0) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(4, 6, 4, 1)
  )

pA <- (pA_donut | pA_meta | pA_prop | pA_n) +
  plot_layout(widths = c(1.25, 1.05, 0.90, 0.80))

# ==============================================================================
# PANEL B: True Raincloud Benign-to-Malignant Gradient (PRJNA658160, n=935)
# ==============================================================================
s658 <- fread(file.path(hgmt_dir, "PRJNA658160_sample_level.csv"))
s658 <- s658[stage %in% c("Health", "Breast Diseases", "Breast Neoplasms")]
s658[, Stage := factor(stage, levels = c("Health", "Breast Diseases", "Breast Neoplasms"),
                       labels = c("Healthy\n(n=430)", "Benign\n(n=108)", "BC\n(n=397)"))]
grad <- fread(file.path(res_dir, "HGMT_benign_to_malignant_gradient.csv"))[project_id == "PRJNA658160"]

b_df <- rbind(
  data.table(Stage = s658$Stage, feature = "Richness", val = s658$Richness,
             feat_lbl = "'Taxon Richness'"),
  data.table(Stage = s658$Stage, feature = "Shannon", val = s658$Shannon,
             feat_lbl = "'Shannon Diversity'"),
  data.table(Stage = s658$Stage, feature = "Lachnospira", val = log10(s658$Lachnospira + 0.01),
             feat_lbl = "italic('Lachnospira')~'(log'[10]*' %)'"),
  data.table(Stage = s658$Stage, feature = "Roseburia", val = log10(s658$Roseburia + 0.01),
             feat_lbl = "italic('Roseburia')~'(log'[10]*' %)'")
)
feat_order <- c("'Taxon Richness'", "'Shannon Diversity'",
                "italic('Lachnospira')~'(log'[10]*' %)'", "italic('Roseburia')~'(log'[10]*' %)'")
b_df[, feat_f := factor(feat_lbl, levels = feat_order)]
b_df[, stage_x := as.numeric(Stage)]

# Build smooth left half-violin polygons manually per (feat_f, Stage) for true raincloud look
dens_rows <- list()
for (ff in feat_order) {
  sub_ff <- b_df[feat_f == ff]
  rng_y <- range(sub_ff$val, na.rm = TRUE)
  for (st in levels(s658$Stage)) {
    sub_s <- sub_ff[Stage == st]
    sx <- match(st, levels(s658$Stage))
    d <- density(sub_s$val, from = rng_y[1], to = rng_y[2], n = 128, na.rm = TRUE)
    w <- d$y / max(d$y) * 0.36
    # Polygon on left side of sx
    poly_x <- c(rep(sx - 0.04, length(d$x)), rev(sx - 0.04 - w))
    poly_y <- c(d$x, rev(d$x))
    dens_rows[[length(dens_rows) + 1]] <- data.table(
      feat_f = factor(ff, levels = feat_order),
      Stage = factor(st, levels = levels(s658$Stage)),
      x = poly_x, y = poly_y
    )
  }
}
dens_dt <- rbindlist(dens_rows)

set.seed(2026)
b_df[, jx := stage_x + runif(.N, 0.07, 0.28)]
b_means <- b_df[, .(med = median(val, na.rm = TRUE), stage_x = stage_x[1]), by = .(feat_f, Stage)]

stat_ann <- merge(unique(b_df[, .(feature, feat_f)]), grad[, .(feature, spearman_rho, spearman_p)], by = "feature")
yranges <- b_df[, .(ymax = max(val, na.rm = TRUE), ymin = min(val, na.rm = TRUE)), by = feat_f]
stat_ann <- merge(stat_ann, yranges, by = "feat_f")
stat_ann[, lbl := sprintf("rho == %.2f*','~~italic(P) == %s", spearman_rho,
                          ifelse(spearman_p < 1e-3,
                                 sprintf("%.1f%%*%%10^{%d}",
                                         as.numeric(strsplit(formatC(spearman_p, format="e", digits=1), "e")[[1]][1]),
                                         as.integer(strsplit(formatC(spearman_p, format="e", digits=1), "e")[[1]][2])),
                                 sprintf("%.3f", spearman_p))), by = feature]

pB <- ggplot() +
  geom_polygon(data = dens_dt, aes(x = x, y = y, fill = Stage, color = Stage), alpha = 0.28, linewidth = 0.35) +
  geom_point(data = b_df, aes(x = jx, y = val, color = Stage), size = 0.38, alpha = 0.28) +
  geom_boxplot(data = b_df, aes(x = stage_x - 0.04, y = val, fill = Stage, group = Stage),
               width = 0.13, outlier.shape = NA, alpha = 0.90, color = PAL$ink, linewidth = 0.35) +
  geom_line(data = b_means, aes(x = stage_x - 0.04, y = med, group = 1),
            color = PAL$ink, linewidth = 0.7) +
  geom_point(data = b_means, aes(x = stage_x - 0.04, y = med),
             color = "white", fill = PAL$ink, shape = 21, size = 2.0, stroke = 0.5) +
  geom_text(data = stat_ann, aes(x = 3.38, y = ymax + 0.16 * (ymax - ymin), label = lbl),
            parse = TRUE, hjust = 1, vjust = 1, size = 2.35, family = "Arial", color = PAL$ink) +
  facet_wrap(~feat_f, scales = "free_y", ncol = 2, labeller = label_parsed) +
  scale_fill_manual(values = c(PAL$hc, PAL$benign, PAL$bc), guide = "none") +
  scale_color_manual(values = c(PAL$hc, PAL$benign, PAL$bc), guide = "none") +
  scale_x_continuous(breaks = 1:3, labels = levels(s658$Stage)) +
  scale_y_continuous(expand = expansion(mult = c(0.06, 0.22))) +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    strip.text = element_text(face = "bold", size = 7.9, margin = margin(b = 4, t = 4)),
    axis.text.x = element_text(size = 7.0),
    plot.margin = margin(14, 6, 4, 6)
  )

# ==============================================================================
# PANEL C: Taxonomic Family Strip + Multi-Cohort Matrix + Range-Whisker Track
# ==============================================================================
g5 <- fread(file.path(res_dir, "HGMT_5cohorts_genus_diff.csv"))
gm <- fread(file.path(res_dir, "GMrepo_genus_differential_by_cohort.csv"))

fam_map <- c(
  "Lachnospira"    = "Lachnospiraceae",
  "Roseburia"      = "Lachnospiraceae",
  "Merdisoma"      = "Lachnospiraceae",
  "UMGS1375"       = "Lachnospiraceae",
  "Howardella"     = "Lachnospiraceae",
  "Enterocloster"  = "Lachnospiraceae",
  "Ruminococcus_E" = "Acutalibacteraceae",
  "CAG-177"        = "Acutalibacteraceae",
  "CAG-269"        = "CAG-508",
  "Collinsella"    = "Actinobacteriota",
  "Eggerthella"    = "Actinobacteriota",
  "Parasutterella" = "Proteobacteria",
  "Escherichia"    = "Proteobacteria",
  "Turicibacter"   = "Bacilli",
  "Streptococcus"  = "Bacilli"
)
fam_cols <- c(
  "Lachnospiraceae"    = "#2B7A9E",
  "Acutalibacteraceae" = "#1F8A70",
  "CAG-508"            = "#7FB3D5",
  "Actinobacteriota"   = "#D9822B",
  "Proteobacteria"     = "#C84B31",
  "Bacilli"            = "#8E6C8A"
)

sel_genera <- names(fam_map)
gm_agg <- gm[genus %in% sel_genera, .(
  project_id = "GMrepo\n(Meta)",
  log2FC = mean(log2FC, na.rm = TRUE),
  wilcox_p = min(c(wilcox_p, linda_p), na.rm = TRUE),
  wilcox_fdr = min(c(wilcox_fdr, linda_fdr), na.rm = TRUE)
), by = genus]

hg_sub <- g5[genus %in% sel_genera, .(genus, project_id, log2FC, wilcox_p, wilcox_fdr)]
mat_c <- rbind(hg_sub, gm_agg, fill = TRUE)

ord_c <- g5[genus %in% sel_genera, .(
  m_fc = mean(log2FC, na.rm = TRUE),
  med_fc = median(log2FC, na.rm = TRUE),
  min_fc = min(log2FC, na.rm = TRUE),
  max_fc = max(log2FC, na.rm = TRUE)
), by = genus][order(med_fc, m_fc)]
ord_genera <- ord_c$genus

coh_levels <- c("PRJNA658160", "PRJNA383849", "PRJNA453965", "PRJNA345373", "PRJNA751260", "GMrepo\n(Meta)")
coh_labels <- c("PRJNA658160\n(16S, Ghana)", "PRJNA383849\n(16S)", "PRJNA453965\n(WGS, China)",
                "PRJNA345373\n(16S)", "PRJNA751260\n(16S, China)", "GMrepo\n(3 Cohorts)")

grid_c <- CJ(genus = ord_genera, project_id = coh_levels)
mat_c <- merge(grid_c, mat_c, by = c("genus", "project_id"), all.x = TRUE)
mat_c[, genus_f := factor(genus, levels = ord_genera)]
mat_c[, coh_f := factor(project_id, levels = coh_levels, labels = coh_labels)]
mat_c[, fc_clip := pmax(pmin(log2FC, 2.5), -2.5)]
mat_c[, neglogP := pmin(-log10(pmax(wilcox_p, 1e-12)), 8)]
mat_c[, sig_tier := factor(
  fifelse(is.na(wilcox_p), "Absent",
  fifelse(wilcox_fdr < 0.05, "FDR < 0.05",
  fifelse(wilcox_p < 0.05, "P < 0.05", "NS"))),
  levels = c("FDR < 0.05", "P < 0.05", "NS", "Absent")
)]

fam_strip_dt <- data.table(
  genus_f = factor(ord_genera, levels = ord_genera),
  Family  = factor(fam_map[ord_genera], levels = names(fam_cols))
)
italic_tax_labels <- parse(text = sprintf("italic('%s')", gsub("_", " ", ord_genera)))

pC_fam <- ggplot(fam_strip_dt, aes(x = 1, y = genus_f, fill = Family)) +
  geom_tile(color = "white", linewidth = 0.5, width = 0.85) +
  scale_fill_manual(values = fam_cols, name = "Taxon Clade") +
  scale_y_discrete(labels = italic_tax_labels) +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7.6, color = PAL$ink),
    legend.position = "bottom",
    legend.key.size = unit(6.5, "pt"),
    plot.margin = margin(8, 0, 4, 4)
  ) +
  guides(fill = guide_legend(nrow = 2, title.position = "top"))

pC_mat <- ggplot(mat_c, aes(x = coh_f, y = genus_f)) +
  geom_tile(fill = "#F8FAFC", color = "#E2E8F0", linewidth = 0.4) +
  geom_vline(xintercept = 5.5, linetype = "dashed", color = "#64748B", linewidth = 0.55) +
  geom_text(data = mat_c[sig_tier == "Absent"], label = "\u2014", color = "#CBD5E0", size = 2.7) +
  geom_point(data = mat_c[sig_tier != "Absent"],
             aes(fill = fc_clip, size = neglogP, color = sig_tier), shape = 21, stroke = 0.68) +
  scale_fill_gradient2(low = PAL$div_lo, mid = PAL$div_mi, high = PAL$div_hi, midpoint = 0,
                       limits = c(-2.5, 2.5), name = expression(bold(log[2]*"FC"))) +
  scale_color_manual(values = c("FDR < 0.05" = PAL$ink, "P < 0.05" = "#4A5568", "NS" = "#CBD5E0"),
                     drop = TRUE, name = "Significance") +
  scale_size_continuous(range = c(1.6, 4.5), breaks = c(1.3, 3, 6),
                        labels = c("0.05", expression(10^-3), expression(10^-6)),
                        name = expression(bold(italic(P)~"value"))) +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.y = element_blank(),
    axis.text.x = element_text(angle = 28, hjust = 1, size = 6.9, color = PAL$ink),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(7.5, "pt"),
    plot.margin = margin(8, 2, 4, 1)
  ) +
  guides(
    fill = guide_colorbar(barwidth = unit(48, "pt"), barheight = unit(5.5, "pt"), title.position = "top"),
    size = guide_legend(title.position = "top", nrow = 1),
    color = guide_legend(title.position = "top", nrow = 1, override.aes = list(fill = "#94A3B8", size = 3.0))
  )

# Right track: Cross-cohort Min-Max Whisker + Median Diamond + Bar
ord_c[, genus_f := factor(genus, levels = ord_genera)]
ord_c[, Dir := ifelse(med_fc < 0, "Depleted in BC", "Enriched in BC")]

pC_range <- ggplot(ord_c, aes(y = genus_f, color = Dir)) +
  geom_vline(xintercept = 0, linetype = "dashed", color = "#718096", linewidth = 0.4) +
  geom_segment(aes(x = min_fc, xend = max_fc, yend = genus_f), linewidth = 0.7, alpha = 0.45) +
  geom_segment(aes(x = 0, xend = med_fc, yend = genus_f), linewidth = 1.05) +
  geom_point(aes(x = med_fc, fill = Dir), shape = 23, size = 2.5, color = "white", stroke = 0.4) +
  scale_color_manual(values = c("Depleted in BC" = PAL$div_lo, "Enriched in BC" = PAL$div_hi), guide = "none") +
  scale_fill_manual(values = c("Depleted in BC" = PAL$div_lo, "Enriched in BC" = PAL$div_hi), guide = "none") +
  scale_x_continuous(limits = c(-8.2, 3.6), breaks = c(-6, -3, 0, 3)) +
  labs(x = expression(bold("Median & Range "*log[2]*"FC")), y = NULL) +
  theme_top(8.1) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    plot.margin = margin(8, 6, 4, 2)
  )

pC <- (pC_fam | pC_mat | pC_range) +
  plot_layout(widths = c(0.14, 1.85, 0.95), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(6, "pt"),
        legend.margin = margin(t = 0, b = 0))

# ==============================================================================
# PANEL D: Menopause-Stratified WGS Taxa & KEGG Pathways (PRJNA453965)
# ==============================================================================
meno <- fread(file.path(res_dir, "HGMT_PRJNA453965_WGS_menopause_stratified.csv"))
pw_labels <- c(
  "Eggerthella"   = "italic('Eggerthella')",
  "Escherichia"   = "italic('Escherichia')",
  "Turicibacter"  = "italic('Turicibacter')",
  "Roseburia"     = "italic('Roseburia')",
  "Lachnospira"   = "italic('Lachnospira')",
  "map00984"      = "'Steroid degradation (map00984)'",
  "map01040"      = "'Unsat. FA biosynth. (map01040)'",
  "map00380"      = "'Tryptophan metab. (map00380)'",
  "map00071"      = "'Fatty acid degrad. (map00071)'",
  "map00650"      = "'Butanoate metab. (map00650)'",
  "map00220"      = "'Arginine biosynth. (map00220)'"
)
msub <- meno[feature %in% names(pw_labels)]
msub[, stars := fifelse(wilcox_p < 0.001, "***", fifelse(wilcox_p < 0.01, "**", fifelse(wilcox_p < 0.05, "*", "")))]
msub[, Stratum := ifelse(stratum == "Postmenopause", "Postmenopause\n(43 BC vs 46 HC)", "Premenopause\n(19 BC vs 25 HC)")]
msub[, Type := factor(ifelse(type == "Genus", "Bacterial Genus", "KEGG Pathway"),
                      levels = c("Bacterial Genus", "KEGG Pathway"))]

ord_d <- msub[stratum == "Postmenopause"][order(log2FC)]$feature
msub[, feat_f := factor(feature, levels = ord_d)]
msub[, txt_col := ifelse(abs(log2FC) > 1.35, "white", PAL$ink)]

pD <- ggplot(msub, aes(x = Stratum, y = feat_f, fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.75) +
  geom_text(aes(label = sprintf("%+.2f%s", log2FC, ifelse(stars != "", paste0(" ", stars), "")),
                color = txt_col, fontface = ifelse(stars != "", "bold", "plain")),
            size = 2.5, family = "Arial") +
  facet_grid(Type ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_fill_gradient2(low = PAL$div_lo, mid = PAL$div_mi, high = PAL$div_hi, midpoint = 0,
                       limits = c(-2.0, 3.1), name = expression(bold(log[2]*"FC"))) +
  scale_color_identity() +
  scale_y_discrete(labels = function(x) parse(text = pw_labels[x]), position = "right") +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 7.3),
    axis.text.y.right = element_text(size = 7.3, hjust = 0),
    strip.placement = "outside",
    strip.background = element_rect(fill = "#EDF2F7", color = NA),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 7.3, color = PAL$navy),
    legend.position = "bottom",
    legend.margin = margin(t = -3),
    plot.margin = margin(14, 6, 4, 6)
  ) +
  guides(fill = guide_colorbar(barwidth = unit(105, "pt"), barheight = unit(5.5, "pt"),
                               title.position = "left", title.vjust = 0.85))

# ==============================================================================
# PANEL E: Continuous Signed -log10(P) KEGG Bubble Lollipop + Category Strip
# ==============================================================================
pw_diff <- fread(file.path(res_dir, "HGMT_PRJNA453965_WGS_KEGG_pathway_diff.csv"))
pw_sel <- pw_diff[!is.na(pathway_name) & pathway_name != "" & wilcox_fdr < 0.05 & n_kos >= 5 &
                  !grepl("Virion|Carcinogenesis|Hepatocellular|Alcoholic", pathway_name, ignore.case = TRUE)][order(wilcox_p)]
pw_sel <- unique(rbind(pw_sel[1:12], pw_sel[pathway == "map00380"]), by = "pathway")
pw_sel[, clean_name := gsub(" - .*", "", pathway_name)]
pw_sel[, clean_name := gsub(" pathway of hypothalamus-pituitary axis", " (HPA axis)", clean_name)]
pw_sel[, clean_name := gsub("Metabolism of xenobiotics by cytochrome P450", "Xenobiotics metab. (CYP450)", clean_name)]
pw_sel[, clean_name := gsub("Drug metabolism", "Drug metab. (CYP450)", clean_name)]
pw_sel[, clean_name := gsub("Pinene, camphor and geraniol degradation", "Pinene, camphor & geraniol degrad.", clean_name)]

# Continuous signed -log10(P) eliminates artificial BH-FDR tie stacking!
pw_sel[, signed_p := sign(log2FC) * -log10(wilcox_p)]
pw_sel[, Dir := ifelse(log2FC > 0, "Enriched in BC", "Depleted in BC")]
# Functional class tags
pw_sel[, KClass := fifelse(pathway %in% c("map00984", "map00907"), "Steroid / Terpenoid",
                    fifelse(pathway %in% c("map00980", "map00982", "map00364"), "Xenobiotics",
                    fifelse(pathway %in% c("map00380", "map00220", "map02060"), "Amino Acid / PTS", "Cellular / Other")))]
kclass_cols <- c("Steroid / Terpenoid" = "#C84B31", "Xenobiotics" = "#D9822B",
                 "Amino Acid / PTS" = "#1F8A70", "Cellular / Other" = "#4A7C9B")

pw_ord <- pw_sel[order(signed_p)]
pw_sel[, name_f := factor(clean_name, levels = pw_ord$clean_name)]
fdr_cut <- -log10(max(pw_sel$wilcox_p))

pE_cls <- ggplot(pw_sel, aes(x = 1, y = name_f, fill = KClass)) +
  geom_tile(color = "white", linewidth = 0.5, width = 0.85) +
  scale_fill_manual(values = kclass_cols, name = "Module") +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_blank(),
    axis.text.y = element_text(size = 7.2, color = PAL$ink),
    legend.position = "bottom",
    legend.key.size = unit(6.5, "pt"),
    legend.margin = margin(t = -3),
    plot.margin = margin(10, 0, 4, 4)
  ) +
  guides(fill = guide_legend(nrow = 2, title.position = "top"))

pE_lolli <- ggplot(pw_sel, aes(x = signed_p, y = name_f)) +
  geom_vline(xintercept = 0, color = PAL$line, linewidth = 0.45) +
  geom_vline(xintercept = c(-fdr_cut, fdr_cut), linetype = "dotted", color = "#94A3B8", linewidth = 0.45) +
  geom_segment(aes(x = 0, xend = signed_p, yend = name_f, color = Dir), linewidth = 0.85, alpha = 0.8) +
  geom_point(aes(size = n_kos, fill = Dir), shape = 21, color = PAL$ink, stroke = 0.45) +
  geom_text(aes(label = sprintf("%+.2f", log2FC),
                hjust = ifelse(signed_p >= 0, -0.38, 1.38)),
            size = 2.3, fontface = "bold", family = "Arial", color = PAL$ink) +
  scale_fill_manual(values = c("Enriched in BC" = PAL$bc, "Depleted in BC" = PAL$hc), name = NULL) +
  scale_color_manual(values = c("Enriched in BC" = PAL$bc, "Depleted in BC" = PAL$hc), guide = "none") +
  scale_size_continuous(range = c(2.3, 5.0), breaks = c(5, 20, 50), name = "KOs") +
  scale_x_continuous(limits = c(-4.85, 5.95), breaks = c(-4, -2, 0, 2, 4)) +
  labs(x = expression(bold("Signed Significance [sign(FC) \u00D7 "*'-log'[10]*"("*italic(P)*")]")), y = NULL) +
  theme_top(8.1) +
  theme(
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(7.5, "pt"),
    legend.margin = margin(t = -3),
    plot.margin = margin(10, 6, 4, 1)
  ) +
  guides(fill = guide_legend(nrow = 2, override.aes = list(size = 3.0)), size = guide_legend(nrow = 1))

pE <- (pE_cls | pE_lolli) +
  plot_layout(widths = c(0.12, 1.8), guides = "collect") &
  theme(legend.position = "bottom", legend.box = "horizontal", legend.spacing.x = unit(6, "pt"))

# ==============================================================================
# PANEL F: Harmonized Cross-Cohort WGS Conservation (China / France / USA)
# ==============================================================================
wgs3 <- fread(file.path(res_dir, "HGMT_3WGS_BC_cross_cohort_profile.csv"))
wgs3[, Section := factor(ifelse(feature_type == "Genus", "Core Taxa", "KEGG Pathways"),
                         levels = c("Core Taxa", "KEGG Pathways"))]
feat_map_f <- c(
  "Faecalibacterium" = "italic('Faecalibacterium')",
  "Bacteroides"      = "italic('Bacteroides')",
  "Roseburia"        = "italic('Roseburia')",
  "Lachnospira"      = "italic('Lachnospira')",
  "Escherichia"      = "italic('Escherichia')",
  "Eggerthella"      = "italic('Eggerthella')",
  "Turicibacter"     = "italic('Turicibacter')",
  "map00220"         = "'Arginine biosynth.'",
  "map00650"         = "'Butanoate metab.'",
  "map00071"         = "'Fatty acid degrad.'",
  "map00380"         = "'Tryptophan metab.'",
  "map01040"         = "'Unsat. FA biosynth.'",
  "map00984"         = "'Steroid degrad.'"
)
wgs3 <- wgs3[feature %in% names(feat_map_f)]
wgs3[, feat_f := factor(feature, levels = rev(names(feat_map_f)))]
wgs3[, coh_f := factor(project_id,
                       levels = c("PRJNA453965", "PRJNA718520", "PRJNA494824"),
                       labels = c("China WGS\n(n=62)", "France WGS\n(n=76)", "USA WGS\n(n=8)"))]
wgs3[, log_mean := log10(mean + 1e-4)]
wgs3[, val_lbl := ifelse(mean >= 1, sprintf("%.1f", mean),
                  ifelse(mean >= 0.05, sprintf("%.2f", mean), sprintf("%.3f", mean)))]

pF <- ggplot(wgs3, aes(x = coh_f, y = feat_f)) +
  geom_tile(fill = "#F8FAFC", color = "#E2E8F0", linewidth = 0.35) +
  geom_point(aes(size = prevalence * 100, fill = log_mean), shape = 21, color = PAL$ink, stroke = 0.4) +
  geom_text(aes(label = val_lbl), nudge_x = 0.28, hjust = 0, size = 2.15, family = "Arial", color = PAL$muted) +
  facet_grid(Section ~ ., scales = "free_y", space = "free_y", switch = "y") +
  scale_size_continuous(range = c(1.3, 5.2), limits = c(0, 100), breaks = c(10, 50, 100), name = "Prev (%)") +
  scale_fill_gradientn(colors = c("#EDF2F7", "#A3D2CA", "#2B7A9E", "#1E3A5F"),
                       limits = c(-3.2, 2.0), name = expression(bold(log[10]*"(Mean)"))) +
  scale_y_discrete(labels = function(x) parse(text = feat_map_f[x])) +
  scale_x_discrete(expand = expansion(add = c(0.55, 0.72))) +
  labs(x = NULL, y = NULL) +
  theme_top(8.1) +
  theme(
    axis.line = element_blank(),
    axis.ticks = element_blank(),
    axis.text.x = element_text(face = "bold", size = 7.3),
    axis.text.y = element_text(size = 7.3),
    strip.placement = "outside",
    strip.background = element_rect(fill = "#EDF2F7", color = NA),
    strip.text.y.left = element_text(angle = 90, face = "bold", size = 7.3, color = PAL$navy),
    legend.position = "bottom",
    legend.box = "horizontal",
    legend.key.size = unit(8, "pt"),
    legend.margin = margin(t = -3),
    plot.margin = margin(14, 6, 4, 6)
  ) +
  guides(
    fill = guide_colorbar(barwidth = unit(75, "pt"), barheight = unit(5.5, "pt"), title.position = "top"),
    size = guide_legend(title.position = "top", nrow = 1)
  )

# ==============================================================================
# Clean Composite Assembly
# ==============================================================================
wrap_panel <- function(plot_obj, tag_letter, title_txt) {
  wrap_elements(full = plot_obj) +
    labs(tag = tag_letter, title = title_txt) +
    theme(
      plot.tag = element_text(family = "Arial", face = "bold", size = 14, color = "#111111"),
      plot.tag.position = c(0.01, 0.986),
      plot.title = element_text(family = "Arial", face = "bold", size = 9.8, color = "#1A202C",
                                hjust = 0, margin = margin(l = 22, t = 2, b = 10)),
      plot.margin = margin(6, 5, 4, 5)
    )
}

wA <- wrap_panel(pA, "A", "HGMT Breast Cancer Cohorts (N = 1,629 Runs)")
wB <- wrap_panel(pB, "B", "Benign-to-Malignant Gradient (PRJNA658160)")
wC <- wrap_panel(pC, "C", "Cross-Cohort Replicated Dysbiosis")
wD <- wrap_panel(pD, "D", "Menopause-Stratified WGS (PRJNA453965)")
wE <- wrap_panel(pE, "E", "Altered Microbial KEGG Pathways")
wF <- wrap_panel(pF, "F", "Cross-Cohort WGS Conservation")

fig1_v3 <- (wA | wB) / (wC | wD) / (wE | wF) +
  plot_layout(heights = c(0.98, 1.28, 1.15), widths = c(1.12, 1.0))

out_png <- file.path(fig_dir, "Figure1_v2_pilot.png")
out_pdf <- file.path(fig_dir, "Figure1_v2_pilot.pdf")

agg_png(out_png, width = 15.2, height = 15.0, units = "in", res = 300, background = "white")
print(fig1_v3)
dev.off()

cairo_pdf(out_pdf, width = 15.2, height = 15.0)
print(fig1_v3)
dev.off()

cat("SAVED upgraded Figure1_v2_pilot.png and Figure1_v2_pilot.pdf successfully!\n")
