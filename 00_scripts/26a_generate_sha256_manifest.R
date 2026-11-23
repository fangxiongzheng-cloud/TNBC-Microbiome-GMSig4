#!/usr/bin/env Rscript
# Generate SHA-256 provenance manifest for frozen GMSig-4 model and Stage A artifacts.
# Never inspects or hashes any external validation cohort files.
# Run: Rscript --vanilla 00_scripts/26a_generate_sha256_manifest.R [project_root]
.libPaths(c("/root/miniconda3/envs/bioinfo-r-core/lib/R/library", "/root/R/tnbc_microbiome_lib", .libPaths()))
suppressPackageStartupMessages({library(digest)})

args <- commandArgs(trailingOnly = TRUE)
root <- if (length(args)) args[1] else "/mnt/d/Bioinformatics/TNBC_Microbiome"

scripts_dir <- file.path(root, "00_scripts")
res_dir     <- file.path(root, "02_results")
freeze_dir  <- file.path(res_dir, "JTM_upgrade/model_freeze")
out_csv     <- file.path(freeze_dir, "sha256_manifest.csv")

collect_group <- function(group_name, paths) {
  paths <- sort(unique(paths[file.exists(paths) & !dir.exists(paths)]))
  if (length(paths) == 0L) stop("No files matched for required manifest group: ", group_name)
  data.frame(
    artifact_group = group_name,
    relative_path  = sub(paste0("^", root, "/?"), "", paths),
    bytes          = as.numeric(file.info(paths)$size),
    sha256         = vapply(paths, function(p) digest::digest(file = p, algo = "sha256"), character(1)),
    modified_utc   = format(file.info(paths)$mtime, tz = "UTC", usetz = TRUE),
    stringsAsFactors = FALSE
  )
}

f_model <- file.path(res_dir, "GMSig4_Model_Objects.rds")
f_19a   <- list.files(scripts_dir, pattern = "^19a", full.names = TRUE)
f_19bf  <- list.files(scripts_dir, pattern = "^19[b-f]", full.names = TRUE)
f_ps    <- list.files(res_dir, pattern = "^GMSig4_PatientScores_.*\\.csv$", full.names = TRUE)
f_26a   <- list.files(scripts_dir, pattern = "^26a", full.names = TRUE)
f_mf    <- setdiff(list.files(freeze_dir, full.names = TRUE), out_csv)

manifest <- rbind(
  collect_group("prior_model_object", f_model),
  collect_group("prior_prep_scripts_19a", f_19a),
  collect_group("prior_model_scripts_19b_19f", f_19bf),
  collect_group("prior_patient_scores", f_ps),
  collect_group("stage_a_freeze_scripts_26a", f_26a),
  collect_group("model_freeze_artifacts", f_mf)
)

rownames(manifest) <- NULL
stopifnot(nrow(manifest) > 0L, !anyNA(manifest$sha256), all(grepl("^[0-9a-f]{64}$", manifest$sha256)))

write.csv(manifest, out_csv, row.names = FALSE)
cat("Wrote", out_csv, "with", nrow(manifest), "entries.\n")
print(manifest[, c("artifact_group", "relative_path", "bytes", "sha256")])
