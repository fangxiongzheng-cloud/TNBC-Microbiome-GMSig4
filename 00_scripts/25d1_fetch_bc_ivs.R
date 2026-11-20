project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

library(data.table)
library(TwoSampleMR)
library(ieugwasr)

res_dir <- file.path(project_root, "02_results")
ao <- as.data.table(ieugwasr::gwasinfo())
cat("Check NK4A136 in OpenGWAS:\n")
print(ao[grepl("NK4A136", trait, ignore.case = TRUE), .(id, trait)])

cat("Extracting clumped genome-wide significant IVs (p < 5e-8) for Overall BC (1126), ER- (1128), ER- Onco (1135)...\n")
bc_ivs <- extract_instruments(outcomes = c("ieu-a-1126", "ieu-a-1128", "ieu-a-1135"),
                              p1 = 5e-8, clump = TRUE, r2 = 0.001, kb = 10000)
cat("Total extracted BC IVs:", nrow(bc_ivs), "\n")
print(as.data.table(bc_ivs)[, .N, by = id.exposure])
fwrite(bc_ivs, file.path(res_dir, "ReverseMR_BC_Instruments_OpenGWAS.csv"))
