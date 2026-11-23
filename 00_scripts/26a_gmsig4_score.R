# Frozen GMSig-4 scoring implementation. No model fitting or outcome access.
# Input: numeric matrix with genes in rows, samples in columns, already on
# the cohort's documented expression scale. No automatic log transformation.
# The historical model uses gene-wise cohort mean and sample SD (n-1).
score_gmsig4 <- function(expr, weights, samples = colnames(expr)) {
  genes <- c("IFNG", "IDO1", "CXCL10", "ENO1")
  if (!is.matrix(expr) || !is.numeric(expr)) stop("expr must be a numeric matrix")
  if (is.null(rownames(expr)) || anyDuplicated(rownames(expr))) stop("unique gene row names required; aggregate probes first")
  if (is.null(colnames(expr)) || anyDuplicated(colnames(expr))) stop("unique sample column names required")
  if (is.null(names(weights)) || !setequal(names(weights), genes)) stop("four named frozen coefficients required")
  if (any(!is.finite(weights))) stop("nonfinite coefficients")
  if (!all(genes %in% rownames(expr))) stop("all four genes are required")
  if (length(samples) < 2L || anyDuplicated(samples) || !all(samples %in% colnames(expr))) stop("at least two unique known samples required")
  x <- t(expr[genes, samples, drop = FALSE])
  if (any(!is.finite(x))) stop("missing/nonfinite expression; resolve by prespecified protocol, never impute silently")
  gene_sd <- apply(x, 2L, sd)
  if (any(!is.finite(gene_sd) | gene_sd <= 0)) stop("zero/undefined gene SD")
  z <- scale(x, center = TRUE, scale = TRUE)
  score <- as.numeric(z %*% weights[genes])
  data.frame(sampleID = samples, GMSig4 = score,
             RiskGroup = ifelse(score > median(score), "High", "Low"),
             stringsAsFactors = FALSE)
}
# This function deliberately fails on undefined inputs. On complete, nonconstant
# inputs it reproduces 19f exactly. High/Low is only a cohort-median display.
# Scores are relative rankings, not absolute survival probabilities or a fixed
# single-patient clinical assay. Do not substitute training means across platforms.
