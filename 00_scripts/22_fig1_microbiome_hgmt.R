project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)

suppressPackageStartupMessages({
  library(data.table)
  library(ggplot2)
  library(patchwork)
})

res_dir <- file.path(project_root, "02_results")
fig_dir <- file.path(project_root, "03_figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

theme_pub <- function(base_size = 9) {
  theme_bw(base_size = base_size) +
    theme(
      panel.grid.minor = element_blank(),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      plot.title = element_text(face = "bold", size = base_size + 1, hjust = 0),
      plot.subtitle = element_text(size = base_size - 0.5, color = "grey30", hjust = 0),
      axis.title = element_text(face = "bold", size = base_size),
      axis.text = element_text(color = "grey15", size = base_size - 1),
      legend.background = element_blank(),
      legend.key = element_blank()
    )
}

# --- Panel 1A: All 8 HGMT Breast Cancer Cohorts (n=1,606) Architecture ---
meta8 <- fread(file.path(res_dir, "HGMT_all8_cohorts_metadata_summary.csv"))
meta8[, label := paste0(project_id, " (", assay_type, ifelse(country != "", paste0(", ", country), ""), ")")]
m8_long <- melt(meta8[, .(label, n_runs_bacteria, `Breast Cancer` = n_BC, `Healthy Control` = n_Health, `Benign Breast Disease` = n_BreastDiseases)],
                id.vars = c("label", "n_runs_bacteria"), variable.name = "Phenotype", value.name = "Count")
m8_long <- m8_long[Count > 0]
m8_long[, Phenotype := factor(Phenotype, levels = c("Breast Cancer", "Benign Breast Disease", "Healthy Control"))]

p1a <- ggplot(m8_long, aes(x = reorder(label, n_runs_bacteria), y = Count, fill = Phenotype)) +
  geom_col(width = 0.72, color = "white", linewidth = 0.25) +
  geom_text(data = meta8, aes(x = reorder(label, n_runs_bacteria), y = n_runs_bacteria, label = paste0("n=", n_runs_bacteria)),
            inherit.aes = FALSE, hjust = -0.15, size = 2.6, fontface = "bold") +
  coord_flip() +
  scale_fill_manual(values = c("Breast Cancer" = "#E64B35", "Benign Breast Disease" = "#F39B7F", "Healthy Control" = "#4DBBD5")) +
  scale_y_continuous(limits = c(0, 1080), expand = c(0, 0)) +
  labs(title = "A. All 8 HGMT Breast Cancer Microbiome Cohorts (Total N=1,606)",
       subtitle = "3 WGS cohorts w/ KEGG KO tables + 5 16S cohorts across China, Ghana, France, USA",
       x = NULL, y = "Sequenced Samples (Runs)", fill = NULL) +
  theme_pub(8.5) +
  theme(legend.position = c(0.75, 0.22))

# --- Panel 1B: Stepwise Dysbiosis Gradient (Health -> Benign Breast Disease -> Malignant BC) ---
grad <- fread(file.path(res_dir, "HGMT_benign_to_malignant_gradient.csv"))
g658 <- grad[project_id == "PRJNA658160" & feature %in% c("Lachnospira", "Roseburia", "Parasutterella", "Eggerthella")]
g_long <- melt(g658[, .(feature, spearman_rho, spearman_p,
                        `Healthy\n(n=430)` = mean_Health,
                        `Benign\n(n=108)` = mean_BreastDiseases,
                        `Malignant BC\n(n=397)` = mean_BC)],
               id.vars = c("feature", "spearman_rho", "spearman_p"), variable.name = "Stage", value.name = "MeanRelAbund")
# Normalize to Healthy baseline = 1.0 for clean cross-genus stepwise fold visualization
g_long[, FoldVsHealth := MeanRelAbund / MeanRelAbund[Stage == "Healthy\n(n=430)"], by = feature]
g_long[, facet_lbl := sprintf("%s (rho=%.2f, p=%.1e)", feature, spearman_rho, spearman_p)]

p1b <- ggplot(g_long, aes(x = Stage, y = FoldVsHealth, group = feature, color = spearman_rho < 0)) +
  geom_hline(yintercept = 1.0, linetype = "dashed", color = "grey60") +
  geom_line(linewidth = 0.9) +
  geom_point(size = 2.6) +
  facet_wrap(~facet_lbl, scales = "free_y", ncol = 2) +
  scale_color_manual(values = c("TRUE" = "#00A087", "FALSE" = "#E64B35"),
                     labels = c("TRUE" = "Depleted along progression", "FALSE" = "Enriched along progression"), name = NULL) +
  labs(title = "B. Stepwise Benign-to-Malignant Dysbiosis (PRJNA658160, n=935)",
       subtitle = "Richness: rho=-0.288 (p=2.37e-19); Shannon: rho=-0.132 (p=5.12e-5)",
       x = NULL, y = "Mean Abundance (Fold vs Healthy)") +
  theme_pub(8.5) +
  theme(legend.position = "bottom", strip.text = element_text(size = 7.2, face = "bold"))

# --- Panel 1C: Cross-Cohort Consensus Differential Genera + GMrepo Anchoring ---
con <- fread(file.path(res_dir, "HGMT_genus_cross_cohort_consistency.csv"))
gm  <- fread(file.path(res_dir, "GMrepo_genus_cross_cohort_consistency.csv"))
sel_genera <- c("Lachnospira", "Roseburia", "Turicibacter", "Ruminococcus_E", "Faecalibacterium",
                "Holdemanella", "Agathobacter", "Blautia_A", "Parasutterella", "Eggerthella",
                "Escherichia", "Streptococcus", "Bacteroides")
csub <- con[genus %in% sel_genera | (n_p05 >= 2 & abs(mean_log2FC) > 0.35)][order(mean_log2FC)]
csub <- unique(rbind(csub[mean_log2FC < 0][1:min(9, .N)], csub[mean_log2FC > 0][1:min(7, .N)]), by = "genus")
csub <- csub[!is.na(genus)]
csub[, Direction := ifelse(mean_log2FC < 0, "Depleted in BC", "Enriched in BC")]
csub[, sig_txt := paste0(n_p05, "/", n_cohorts, " p<0.05", ifelse(!is.na(wgs_p) & wgs_p < 0.05, " (WGS*)", ""))]

p1c <- ggplot(csub, aes(x = reorder(genus, mean_log2FC), y = mean_log2FC, fill = Direction)) +
  geom_col(width = 0.72, color = "white") +
  geom_hline(yintercept = 0, color = "grey30", linewidth = 0.4) +
  geom_text(aes(label = sig_txt, hjust = ifelse(mean_log2FC >= 0, -0.06, 1.06)), size = 2.45) +
  coord_flip() +
  scale_fill_manual(values = c("Depleted in BC" = "#3C5488", "Enriched in BC" = "#DC0000"), name = NULL) +
  scale_y_continuous(limits = c(min(csub$mean_log2FC) - 2.4, max(csub$mean_log2FC) + 2.8)) +
  labs(title = "C. Replicated Case-Control Dysbiosis Across HGMT Cohorts",
       subtitle = "SCFA-producing Lachnospiraceae depletion vs pathobiont enrichment",
       x = NULL, y = "Mean log2(BC / Control)") +
  theme_pub(8.5) +
  theme(legend.position = "bottom")

# --- Panel 1D: Menopause-Stratified WGS Taxa & KEGG Pathways (PRJNA453965, n=133) ---
meno <- fread(file.path(res_dir, "HGMT_PRJNA453965_WGS_menopause_stratified.csv"))
pw_labels <- c(
  "Lachnospira" = "Genus: Lachnospira",
  "Roseburia" = "Genus: Roseburia",
  "Escherichia" = "Genus: Escherichia",
  "Turicibacter" = "Genus: Turicibacter",
  "Eggerthella" = "Genus: Eggerthella",
  "map00984" = "KEGG: Steroid degradation (map00984)",
  "map00380" = "KEGG: Tryptophan metab. (map00380)",
  "map00650" = "KEGG: Butanoate metab. (map00650)",
  "map00220" = "KEGG: Arginine biosynth. (map00220)",
  "map00071" = "KEGG: Fatty acid degrad. (map00071)",
  "map01040" = "KEGG: Unsat. FA biosynth. (map01040)"
)
msub <- meno[feature %in% names(pw_labels)]
msub[, feat_lbl := pw_labels[feature]]
msub[, stars := ifelse(wilcox_p < 0.001, "***", ifelse(wilcox_p < 0.01, "**", ifelse(wilcox_p < 0.05, "*", "")))]
msub[, Stratum := ifelse(stratum == "Postmenopause", "Postmenopause (43 BC vs 46 HC)", "Premenopause (19 BC vs 25 HC)")]

p1d <- ggplot(msub, aes(x = Stratum, y = reorder(feat_lbl, ifelse(stratum == "Postmenopause", log2FC, 0)), fill = log2FC)) +
  geom_tile(color = "white", linewidth = 0.6) +
  geom_text(aes(label = sprintf("%+.2f %s", log2FC, stars)), size = 2.6, fontface = "bold") +
  scale_fill_gradient2(low = "#3B4CC0", mid = "white", high = "#B40426", midpoint = 0, limits = c(-2.2, 3.1), name = "log2FC") +
  scale_x_discrete(labels = c("Postmenopause\n(43 BC vs 46 HC)", "Premenopause\n(19 BC vs 25 HC)")) +
  labs(title = "D. Menopause-Stratified WGS Taxa & KEGG Pathways (PRJNA453965)",
       subtitle = "*p<0.05, **p<0.01, ***p<0.001 (Mann-Whitney U)",
       x = NULL, y = NULL) +
  theme_pub(8.5)

# --- Panel 1E: Top WGS Microbial KEGG Pathways & 3-WGS Cross-Cohort Consensus ---
pw_diff <- fread(file.path(res_dir, "HGMT_PRJNA453965_WGS_KEGG_pathway_diff.csv"))
pw_sel  <- pw_diff[!is.na(pathway_name) & pathway_name != "" & wilcox_fdr < 0.05 & n_kos >= 5 &
                   !grepl("Virion|Carcinogenesis|Hepatocellular|Alcoholic", pathway_name, ignore.case = TRUE)][order(wilcox_p)]
pw_sel  <- pw_sel[1:min(13, .N)]
pw_sel[, clean_name := gsub(" - .*", "", pathway_name)]
pw_sel[, clean_name := gsub(" pathway of hypothalamus-pituitary axis", " (HPA axis)", clean_name)]
pw_sel[, clean_name := gsub("Metabolism of xenobiotics by cytochrome P450", "Xenobiotics metabolism by CYP450", clean_name)]
pw_sel[, clean_name := gsub("Pinene, camphor and geraniol degradation", "Pinene, camphor & geraniol degradation", clean_name)]
pw_sel[, Dir := ifelse(log2FC > 0, "Enriched in BC", "Depleted in BC")]

p1e <- ggplot(pw_sel, aes(x = reorder(clean_name, -log10(wilcox_fdr)), y = -log10(wilcox_fdr), fill = Dir)) +
  geom_col(width = 0.72, color = "white") +
  geom_text(aes(label = sprintf("%+.2f", log2FC)), hjust = -0.12, size = 2.45) +
  coord_flip() +
  scale_fill_manual(values = c("Enriched in BC" = "#E64B35", "Depleted in BC" = "#4DBBD5"), name = NULL) +
  scale_y_continuous(limits = c(0, max(-log10(pw_sel$wilcox_fdr)) * 1.25), expand = c(0, 0)) +
  labs(title = "E. Top Altered WGS KEGG Pathways (PRJNA453965, FDR<0.05)",
       subtitle = "Numbers indicate log2(BC / Healthy Control) pathway KO sum (n_KO >= 5)",
       x = NULL, y = "-log10(BH FDR)") +
  theme_pub(8.5) +
  theme(legend.position = "bottom", axis.text.y = element_text(size = 7.4))

# --- Panel 1F: Multi-Cohort WGS Functional & SCFA-Genus Profile across China (PRJNA453965), France (PRJNA718520), USA (PRJNA494824) ---
wgs3 <- fread(file.path(res_dir, "HGMT_3WGS_BC_cross_cohort_profile.csv"))
wgs3_pw <- wgs3[feature_type == "KEGG_Pathway"]
wgs3_pw[, pw_short := c(
  "map00984" = "Steroid degrad.",
  "map00380" = "Tryptophan metab.",
  "map00071" = "Fatty acid degrad.",
  "map00650" = "Butanoate metab.",
  "map00220" = "Arginine biosynth.",
  "map01040" = "Unsat. FA biosynth."
)[feature]]
wgs3_pw[, coh_lbl := c(
  "PRJNA453965" = "China WGS\n(n=62 BC)",
  "PRJNA718520" = "France WGS\n(n=76 BC)",
  "PRJNA494824" = "USA WGS\n(n=8 BC)"
)[project_id]]

p1f <- ggplot(wgs3_pw, aes(x = coh_lbl, y = pw_short, size = prevalence * 100, color = log10(mean + 1e-4))) +
  geom_point() +
  geom_text(aes(label = sprintf("%.3f", mean)), size = 2.45, color = "grey15", nudge_y = 0.28, vjust = 0) +
  scale_size_continuous(range = c(4, 7.5), name = "Prevalence (%)") +
  scale_color_viridis_c(option = "plasma", name = "log10(Mean KO)") +
  scale_y_discrete(expand = expansion(add = c(0.45, 0.72))) +
  labs(title = "F. Cross-Cohort WGS KEGG Functional Conservation",
       subtitle = "100% prevalence across Chinese, French & US BC microbiomes",
       x = NULL, y = NULL) +
  theme_pub(8.5)

fig1 <- (p1a | p1b) / (p1c | p1d) / (p1e | p1f) +
  plot_annotation(
    title = "Figure 1. Cross-Cohort Breast Cancer Microbiome Dysbiosis, Benign-to-Malignant Gradient, and Functional Pathways (8 HGMT Cohorts, N=1,606)",
    theme = theme(plot.title = element_text(face = "bold", size = 12, hjust = 0))
  )

ggsave(file.path(fig_dir, "Figure1.pdf"), fig1, width = 15.5, height = 14.5, device = cairo_pdf)
ggsave(file.path(fig_dir, "Figure1.png"), fig1, width = 15.5, height = 14.5, dpi = 300)
cat("SAVED Figure1.pdf and Figure1.png successfully!\n")
