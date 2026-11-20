project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

library(data.table)
library(ieugwasr)
library(jsonlite)

res_dir <- file.path(project_root, "02_results")
ao <- as.data.table(ieugwasr::gwasinfo())

# Find all Kurilshikov (MiBioGen) IDs in OpenGWAS matching our replicated causal taxa
mib_kur <- ao[grepl("Kurilshikov", author, ignore.case = TRUE)]
cat("Kurilshikov MiBioGen studies in OpenGWAS:", nrow(mib_kur), "\n")
target_patterns <- c("Dorea", "Sellimonas", "FamilyXIIIUCG001", "Family XIII UCG001",
                     "Desulfovibrio", "LachnospiraceaeNK4A136", "Eubacteriumxylanophilum",
                     "Eubacterium xylanophilum", "Parabacteroides", "Marvinbryantia",
                     "Oxalobacter", "Anaerostipes", "Ruminococcaceae", "Lactobacillales",
                     "Adlercreutzia", "Ruminiclostridium6")
matched <- mib_kur[grepl(paste(target_patterns, collapse = "|"), trait, ignore.case = TRUE)]
print(matched[, .(id, trait, sample_size)])
fwrite(matched[, .(id, trait, sample_size)], file.path(res_dir, "OpenGWAS_MiBioGen_Target_IDs.csv"))

# Search Breast Cancer studies in OpenGWAS (including Zhang / TNBC if present)
bc_all <- ao[grepl("breast cancer", trait, ignore.case = TRUE)]
print(bc_all[grepl("triple|ER-|TNBC|BCAC|Zhang|Michailidou", paste(trait, author, consortium), ignore.case = TRUE),
             .(id, trait, ncase, ncontrol, author)][1:min(30, .N)])

# Query EBI GWAS Catalog REST API for GCST010100 & GCST90454344 genome-wide lead associations if available
for (st in c("GCST010100", "GCST90454344")) {
  url <- sprintf("https://www.ebi.ac.uk/gwas/rest/api/studies/%s/associations?projection=associationByStudy", st)
  js <- tryCatch(fromJSON(url), error = function(e) NULL)
  if (!is.null(js) && !is.null(js$`_embedded`$associations)) {
    cat(st, "GWAS Catalog curated associations:", nrow(js$`_embedded`$associations), "\n")
  } else {
    cat(st, "GWAS Catalog REST API returned 0 or unavailable\n")
  }
}
