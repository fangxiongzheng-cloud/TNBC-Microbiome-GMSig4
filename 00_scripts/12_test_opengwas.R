project_root <- normalizePath(Sys.getenv("PROJECT_ROOT", unset = "."), winslash = "/", mustWork = TRUE)
project_env <- file.path(project_root, ".Renviron")
if (!nzchar(Sys.getenv("OPENGWAS_JWT")) && file.exists(project_env)) readRenviron(project_env)

cat("OPENGWAS_JWT set:", nzchar(Sys.getenv("OPENGWAS_JWT")), "\n")
invisible(ieugwasr::user())
cat("OpenGWAS authentication succeeded.\n")

ids <- c("ieu-a-1126", "ieu-a-1128", "ieu-a-1135", "ebi-a-GCST010098", "ebi-a-GCST010100", "ebi-a-GCST90454344", "ebi-a-GCST90016908")
cat("\n=== Query OpenGWAS gwasinfo ===\n")
info <- tryCatch(ieugwasr::gwasinfo(ids), error = function(e) { cat("err:", e$message, "\n"); NULL })
if (!is.null(info)) {
  print(info[, intersect(c("id", "trait", "ncase", "ncontrol", "sample_size", "nsnp", "pmid", "year"), colnames(info))])
}

# Check also if eQTLGen blood eQTL is on OpenGWAS (eqtl-a-ENSG...)
cat("\n=== Test OpenGWAS eQTLGen & LD clump ===\n")
eq_test <- tryCatch(ieugwasr::gwasinfo("eqtl-a-ENSG00000169245"), error = function(e) { cat("eqtl err:", e$message, "\n"); NULL }) # CXCL10
if (!is.null(eq_test)) print(eq_test[, intersect(c("id", "trait", "sample_size", "nsnp"), colnames(eq_test))])

cl_test <- tryCatch({
  ieugwasr::ld_clump(data.frame(rsid = c("rs10830289", "rs10830290"), pval = c(1e-6, 1e-5), id = "test"))
}, error = function(e) { cat("clump err:", e$message, "\n"); NULL })
if (!is.null(cl_test)) print(cl_test)
