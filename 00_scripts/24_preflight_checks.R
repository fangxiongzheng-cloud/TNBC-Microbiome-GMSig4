project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

library(data.table)
library(pROC)
library(TwoSampleMR)
library(ieugwasr)

res_dir   <- file.path(project_root, "02_results")
cache_dir <- file.path(project_root, "01_data", "TCGA_GEO", "clean_cache")

# 1. Exact weights in GMSig4_Model_Objects.rds + verify GSE194040 AUC
mod <- readRDS(file.path(res_dir, "GMSig4_Model_Objects.rds"))
cat("weights:\n"); print(mod$weights, digits = 10)
cat("uni_betas:\n"); print(mod$uni_betas, digits = 10)

g194 <- readRDS(file.path(cache_dir, "GSE194040_clean.rds"))
g173 <- readRDS(file.path(cache_dir, "GSE173839_clean.rds"))
genes4 <- names(mod$weights)

p_all <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !is.na(g194$clin$pCR))
p_tn  <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
p_hrp <- which(grepl("Pembrolizumab", g194$clin$arm, ignore.case = TRUE) & !g194$clin$is_TNBC & !is.na(g194$clin$pCR))
c_tn  <- which(g194$clin$arm == "Paclitaxel" & g194$clin$is_TNBC & !is.na(g194$clin$pCR))
hr_col <- grep("^HR\\.status", colnames(g173$clin), value = TRUE)[1]
d_tn  <- which(g173$clin$Arm == "durvalumab/olaparib" & g173$clin$pCR.status %in% c(0, 1) & g173$clin[[hr_col]] == 0)

export_subgroup <- function(expr, idx, y01, label) {
  M <- scale(t(expr[genes4, idx, drop = FALSE]))
  dt <- data.table(
    subgroup = label,
    pCR = y01,
    Ridge_Score = as.numeric(M %*% mod$weights),
    UniBeta_Score = as.numeric(M %*% mod$uni_betas)
  )
  cat(sprintf("%-22s n=%d Ridge_AUC=%.6f UniBeta_AUC=%.6f\n",
      label, nrow(dt),
      as.numeric(auc(roc(dt$pCR, -dt$Ridge_Score, quiet = TRUE, direction = "<"))),
      as.numeric(auc(roc(dt$pCR, -dt$UniBeta_Score, quiet = TRUE, direction = "<")))))
  dt
}

sub_scores <- rbindlist(list(
  export_subgroup(g194$expr, p_tn,  g194$clin$pCR[p_tn],  "Pembro_TNBC"),
  export_subgroup(g194$expr, p_all, g194$clin$pCR[p_all], "Pembro_All"),
  export_subgroup(g194$expr, p_hrp, g194$clin$pCR[p_hrp], "Pembro_HRp"),
  export_subgroup(g194$expr, c_tn,  g194$clin$pCR[c_tn],  "Paclitaxel_Ctrl_TNBC"),
  export_subgroup(g173$expr, d_tn,  g173$clin$pCR.status[d_tn], "DurvaOlap_TNBC")
))
fwrite(sub_scores, file.path(res_dir, "GMSig4_ICI_Subgroup_ExactScores.csv"))

# 2. Check package availability for CellChat / dcurves / rmda
for (pkg in c("CellChat", "dcurves", "rmda", "survival")) {
  cat(sprintf("Pkg %-10s installed: %s\n", pkg, requireNamespace(pkg, quietly = TRUE)))
}

# 3. Check OpenGWAS user + search MiBioGen / Kurilshikov / gut microbiome studies
u <- tryCatch(ieugwasr::user(), error = function(e) e)
cat("OpenGWAS user status:", if (inherits(u, "error")) u$message else "OK", "\n")
ao <- tryCatch(as.data.table(ieugwasr::gwasinfo()), error = function(e) data.table())
cat("gwasinfo returned rows:", nrow(ao), "\n")
if (nrow(ao) > 0) {
  mib_ao <- ao[grepl("Dorea|Sellimonas|Desulfovibrio|Marvinbryantia|Anaerostipes|MiBioGen|Kurilshikov", paste(trait, author, consortium), ignore.case = TRUE)]
  cat("Matched microbiome GWAS in OpenGWAS:", nrow(mib_ao), "\n")
  if (nrow(mib_ao) > 0) print(head(mib_ao[, .(id, trait, ncase, ncontrol, sample_size, author)], 20))
  bc_ao <- ao[id %in% c("ieu-a-1126", "ieu-a-1128", "ieu-a-1135") | grepl("triple negative|breast cancer.*Zhang", paste(trait, author), ignore.case = TRUE)]
  cat("Matched BC/TNBC GWAS in OpenGWAS:\n")
  print(head(bc_ao[, .(id, trait, ncase, ncontrol, sample_size, author)], 15))
}
