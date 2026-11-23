#!/usr/bin/env Rscript
# Reconstruct and freeze existing GMSig-4; never fits on external outcomes.
# Run: Rscript --vanilla 00_scripts/26a_freeze_gmsig4.R [project_root]
.libPaths(c("/root/miniconda3/envs/bioinfo-r-core/lib/R/library", "/root/R/tnbc_microbiome_lib", .libPaths()))
suppressPackageStartupMessages({library(survival); library(glmnet); library(jsonlite)})
args <- commandArgs(trailingOnly=TRUE)
root <- if(length(args)) args[1] else "/mnt/d/Bioinformatics/TNBC_Microbiome"
out <- file.path(root,"02_results/JTM_upgrade/model_freeze")
dir.create(out,recursive=TRUE,showWarnings=FALSE)
source(file.path(root,"00_scripts/26a_gmsig4_score.R"))
options(digits=17)
res <- file.path(root,"02_results")
cache <- file.path(root,"01_data/TCGA_GEO/clean_cache")
m <- readRDS(file.path(res,"GMSig4_Model_Objects.rds"))
w <- m$weights
stopifnot(identical(names(w),c("IFNG","IDO1","CXCL10","ENO1")))
write_json(list(model="GMSig-4", authority="Exact double-precision weights from pre-existing GMSig4_Model_Objects.rds; never replaced by reconstruction", genes=names(w), weights=as.list(w), expression_orientation="genes_by_samples", preprocessing="already transformed gene expression; gene-wise cohort centering and sample SD (n-1), then weighted sum", risk_direction="higher score = higher OS hazard", threshold="strictly greater than cohort median = High; otherwise Low; descriptive only", single_sample_validated=FALSE, absolute_risk_validated=FALSE),file.path(out,"frozen_coefficients.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
write.csv(data.frame(gene=names(w),coefficient=as.numeric(w)),file.path(out,"frozen_coefficients.csv"),row.names=FALSE)
cat("Exact saved coefficients exported\n");print(w)
tc <- readRDS(file.path(cache,"TCGA_BRCA_clean.rds"))
idx <- which(tc$clin$is_TNBC & !is.na(tc$clin$OS) & !is.na(tc$clin$OS.time) & tc$clin$OS.time>0)
x <- t(tc$expr[names(w),idx,drop=FALSE]); z <- scale(x)
y <- Surv(tc$clin$OS.time[idx]/365.25,tc$clin$OS[idx])
write.csv(data.frame(gene=names(w),training_mean=colMeans(x),training_sample_sd=apply(x,2,sd)),file.path(out,"training_scaling.csv"),row.names=FALSE)
write.csv(data.frame(sampleID=colnames(tc$expr)[idx],training_order=seq_along(idx)),file.path(out,"training_sample_order.csv"),row.names=FALSE)
saved_scores <- read.csv(file.path(res,"GMSig4_PatientScores_TCGA_TNBC.csv"),check.names=FALSE)
sc <- score_gmsig4(tc$expr,w,colnames(tc$expr)[idx])
mt <- match(sc$sampleID,saved_scores$sampleID)
stopifnot(!anyNA(mt))
tests <- list()
add <- function(name,value,tolerance,passed) tests[[length(tests)+1L]] <<- data.frame(check=name,value=as.character(value),tolerance=as.character(tolerance),passed=passed)
add("TCGA score exact reconstruction max abs error",max(abs(sc$GMSig4-saved_scores$GMSig4[mt])),1e-12,max(abs(sc$GMSig4-saved_scores$GMSig4[mt]))<1e-12)
add("TCGA z-expression max abs error",max(abs(z-as.matrix(saved_scores[mt,names(w)]))),1e-12,max(abs(z-as.matrix(saved_scores[mt,names(w)])))<1e-12)
add("TCGA risk groups identical",all(sc$RiskGroup==saved_scores$RiskGroup[mt]),TRUE,all(sc$RiskGroup==saved_scores$RiskGroup[mt]))
set.seed(2026)
cv <- cv.glmnet(z,y,family="cox",alpha=0,nfolds=5)
wr <- as.numeric(coef(cv,s="lambda.min")); names(wr)<-names(w)
err<-max(abs(w-wr))
add("Ridge reconstruction coefficient max abs error",err,1e-12,err<1e-12)
write.csv(data.frame(gene=names(w),saved=as.numeric(w),reconstructed=wr,difference=wr-w),file.path(out,"ridge_reconstruction_check.csv"),row.names=FALSE)
# Original cv.glmnet makes the same balanced random fold allocation internally.
# Recover it explicitly and verify that specifying foldid reproduces coefficients.
set.seed(2026); foldid<-sample(rep(seq(5),length=nrow(z)))
cv_fold <- cv.glmnet(z,y,family="cox",alpha=0,nfolds=5,foldid=foldid)
add("Explicit fold allocation coefficient match",max(abs(as.numeric(coef(cv_fold,s="lambda.min"))-wr)),1e-12,max(abs(as.numeric(coef(cv_fold,s="lambda.min"))-wr))<1e-12)
write.csv(data.frame(sampleID=rownames(z),fold=foldid),file.path(out,"reconstructed_training_cv_folds.csv"),row.names=FALSE)
saveRDS(cv,file.path(out,"reconstructed_training_cv_object.rds"))
# Archive existing nomogram baseline only; it is NOT the ridge baseline.
# Original objects omit model frame; restore their named data from existing exports,
# without refitting or reading any new validation outcome.
tcga_df <- saved_scores
g960_df <- read.csv(file.path(res,"GMSig4_PatientScores_GSE96058_TNBC_Ext.csv"),check.names=FALSE)
tcga_cc <- tcga_df[complete.cases(tcga_df[,c("OS_yr","OS","GMSig4","Age","Stage_III_IV","Node_Pos")]),]
g960_cc <- g960_df[complete.cases(g960_df[,c("OS_yr","OS","GMSig4","Age","Tumor_mm","Node_Pos")]),]
for(nm in c("nomo_tcga","nomo_g960")) {
  f<-m[[nm]]
  write.csv(basehaz(f,centered=FALSE),file.path(out,paste0(nm,"_existing_baseline_uncentered.csv")),row.names=FALSE)
  write.csv(data.frame(term=names(coef(f)),coefficient=as.numeric(coef(f)),cox_center=f$means),file.path(out,paste0(nm,"_existing_terms.csv")),row.names=FALSE)
}
# Verify scoring under reordered genes/weights; exercise specified fail-fast checks.
reord<-score_gmsig4(tc$expr[rev(names(w)),idx,drop=FALSE],w[rev(names(w))])
add("Gene/weight ordering invariant",max(abs(reord$GMSig4-sc$GMSig4)),1e-12,max(abs(reord$GMSig4-sc$GMSig4))<1e-12)
reject<-function(a) inherits(tryCatch(a,error=function(e)e),"error")
small<-tc$expr[names(w),idx[1:4],drop=FALSE]
add("Missing gene rejected",reject(score_gmsig4(small[-1,,drop=FALSE],w)),TRUE,reject(score_gmsig4(small[-1,,drop=FALSE],w)))
bad<-small;bad[1,1]<-NA_real_
add("Nonfinite input rejected",reject(score_gmsig4(bad,w)),TRUE,reject(score_gmsig4(bad,w)))
bad<-small;bad[1,]<-1
add("Zero gene SD rejected",reject(score_gmsig4(bad,w)),TRUE,reject(score_gmsig4(bad,w)))
add("Single sample rejected",reject(score_gmsig4(small[,1,drop=FALSE],w)),TRUE,reject(score_gmsig4(small[,1,drop=FALSE],w)))
# Other existing patient tables contain already transformed z genes: no new outcomes.
for(nm in c("GSE58812","GSE96058_TNBC_Ext")) {
  p<-read.csv(file.path(res,paste0("GMSig4_PatientScores_",nm,".csv")),check.names=FALSE)
  e<-max(abs(as.numeric(as.matrix(p[,names(w)]) %*% w)-p$GMSig4))
  add(paste(nm,"saved z-expression weighted sum"),e,1e-12,e<1e-12)
}
write.csv(do.call(rbind,tests),file.path(out,"verification_checks.csv"),row.names=FALSE)
manifest <- list(freeze_utc=format(Sys.time(),tz="UTC",usetz=TRUE),model="GMSig-4",source_model="02_results/GMSig4_Model_Objects.rds",source_pipeline="00_scripts/19f_signature_clinical_nomogram.R",weights=as.list(w),training=list(dataset="TCGA-BRCA",n=length(idx),OS_events=sum(tc$clin$OS[idx]),eligibility="Primary tumor, ER/PR/HER2 Negative as resolved in 19a; nonmissing OS/OS.time; OS.time>0",time_unit="years = OS.time / 365.25",gene_means=as.list(colMeans(x)),gene_sample_sd=as.list(apply(x,2,sd)),score_sd=sd(sc$GMSig4),score_median=median(sc$GMSig4)),fit=list(family="cox",penalty="ridge",alpha=0,nfolds=5,seed=2026,lambda_rule="lambda.min",lambda_min_reconstructed=cv$lambda.min,lambda_1se_reconstructed=cv$lambda.1se,coefficient_max_abs_error=err,coefficient_reconstruction_matches=err<1e-12,glmnet_internal_standardize=TRUE,cox_ties="Breslow in glmnet",folds_originally_saved=FALSE,cv_object_originally_saved=FALSE),scoring=list(formula="sum_g beta_g * (x_g - cohort_mean_g) / cohort_sample_sd_g",no_score_rescaling=TRUE,median_split="score > cohort median",missing_or_zero_sd="historical scale() had no explicit safety rule; new wrapper fails and requires prespecified handling",platform_portability="cohort-relative; fixed-reference single-patient use not validated"),baseline=list(original_ridge_baseline_saved=FALSE,original_ridge_baseline_reconstruction="not estimated; would require explicitly chosen estimator/ties handling and validation",existing_nomograms="nomo_tcga and nomo_g960 are separate unpenalized Cox models of GMSig4 plus clinical variables; their baselines exported only for provenance; cannot substitute for ridge baseline"),reporting_issue=list(file="GMSig4_Prognostic_Validation_Summary.csv",column="Continuous_HR_perSD",actual_estimand="HR per one raw weighted-score unit; 19f uses coxph(y ~ risk_score) without scale(risk_score)"),verification_passed=all(vapply(tests,function(d)d$passed,logical(1))),package_versions=as.list(sapply(c("survival","glmnet","Matrix","jsonlite"),function(p)as.character(packageVersion(p)))))
write_json(manifest,file.path(out,"model_manifest.json"),auto_unbox=TRUE,pretty=TRUE,digits=NA)
capture.output(sessionInfo(),file=file.path(out,"R_sessionInfo.txt"))
cat("Training n:",length(idx),"events:",sum(tc$clin$OS[idx]),"lambda.min:",cv$lambda.min,"coefficient max error:",err,"\n")
print(do.call(rbind,tests));stopifnot(manifest$verification_passed)
