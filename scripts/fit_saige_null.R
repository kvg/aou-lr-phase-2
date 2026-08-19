#!/usr/bin/env Rscript
# Fit a SAIGE logistic null model (step 1) for one binary phenotype.
#
# Aligns pheno_cov to analysis_samples, keeps complete cases for the phenotype
# + selected covariates, then calls SAIGE step1_fitNULLGLMM.R with the shared
# sparse GRM / PLINK prefix from build_saige_plink_and_grm.sh.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    pheno_cov = NA_character_,
    phenotype = NA_character_,
    covariates = NA_character_,
    analysis_samples = NA_character_,
    plink_prefix = NA_character_,
    sparse_grm = NA_character_,
    sparse_grm_ids = NA_character_,
    step1_r = "/opt/SAIGE/extdata/step1_fitNULLGLMM.R",
    n_threads = 8L,
    out_prefix = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (i == length(args)) stop(paste("Missing value for", key))
    val <- args[[i + 1]]
    if (key == "--pheno-cov") out$pheno_cov <- val
    else if (key == "--phenotype") out$phenotype <- val
    else if (key == "--covariates") out$covariates <- val
    else if (key == "--analysis-samples") out$analysis_samples <- val
    else if (key == "--plink-prefix") out$plink_prefix <- val
    else if (key == "--sparse-grm") out$sparse_grm <- val
    else if (key == "--sparse-grm-ids") out$sparse_grm_ids <- val
    else if (key == "--step1-r") out$step1_r <- val
    else if (key == "--n-threads") out$n_threads <- as.integer(val)
    else if (key == "--out-prefix") out$out_prefix <- val
    else stop(paste("Unknown arg:", key))
    i <- i + 2
  }
  req <- c("pheno_cov", "phenotype", "covariates", "analysis_samples",
           "plink_prefix", "sparse_grm", "sparse_grm_ids", "out_prefix")
  for (r in req) {
    if (is.na(out[[r]]) || !nzchar(out[[r]])) stop(paste("Missing required", r))
  }
  out
}

opt <- parse_args(args)

covars <- scan(opt$covariates, what = character(), quiet = TRUE)
if (length(covars) < 1) stop("No covariates listed")

samples <- scan(opt$analysis_samples, what = character(), quiet = TRUE)
pheno <- fread(opt$pheno_cov, sep = "\t", header = TRUE, data.table = FALSE)
if (!("ID" %in% names(pheno))) stop("pheno_cov must contain ID column")
if (!(opt$phenotype %in% names(pheno))) {
  stop(sprintf("Phenotype %s not in pheno_cov", opt$phenotype))
}
missing_cov <- setdiff(covars, names(pheno))
if (length(missing_cov) > 0) {
  stop(sprintf("Missing covariate columns: %s", paste(missing_cov, collapse = ",")))
}

pheno$ID <- as.character(pheno$ID)
pheno <- pheno[match(samples, pheno$ID), , drop = FALSE]
if (any(is.na(pheno$ID))) stop("Some analysis_samples missing from pheno_cov")

keep_cols <- c("ID", opt$phenotype, covars)
cc <- stats::complete.cases(pheno[, keep_cols])
pheno <- pheno[cc, , drop = FALSE]
pheno[[opt$phenotype]] <- as.numeric(pheno[[opt$phenotype]])
for (cv in covars) {
  pheno[[cv]] <- as.numeric(pheno[[cv]])
}

# Binary check
vals <- sort(unique(pheno[[opt$phenotype]]))
if (!all(vals %in% c(0, 1))) {
  stop(sprintf("Phenotype %s is not binary 0/1 after filtering: %s",
               opt$phenotype, paste(vals, collapse = ",")))
}

dir.create(dirname(opt$out_prefix), recursive = TRUE, showWarnings = FALSE)
pheno_path <- paste0(opt$out_prefix, ".pheno.tsv")
fwrite(pheno[, keep_cols, drop = FALSE], pheno_path, sep = "\t", na = "NA", quote = FALSE)

sample_path <- paste0(opt$out_prefix, ".samples.txt")
writeLines(pheno$ID, sample_path)

covar_list <- paste(covars, collapse = ",")
meta_path <- paste0(opt$out_prefix, ".null_meta.tsv")
fwrite(
  data.frame(
    phenotype = opt$phenotype,
    n_samples = nrow(pheno),
    n_cases = sum(pheno[[opt$phenotype]] == 1),
    n_controls = sum(pheno[[opt$phenotype]] == 0),
    covariates = covar_list,
    stringsAsFactors = FALSE
  ),
  meta_path,
  sep = "\t",
  quote = FALSE
)

cmd <- paste(
  "Rscript", shQuote(opt$step1_r),
  paste0("--plinkFile=", shQuote(opt$plink_prefix)),
  paste0("--phenoFile=", shQuote(pheno_path)),
  paste0("--phenoCol=", opt$phenotype),
  paste0("--covarColList=", covar_list),
  "--sampleIDColinphenoFile=ID",
  "--traitType=binary",
  paste0("--outputPrefix=", shQuote(opt$out_prefix)),
  paste0("--sparseGRMFile=", shQuote(opt$sparse_grm)),
  paste0("--sparseGRMSampleIDFile=", shQuote(opt$sparse_grm_ids)),
  "--useSparseGRMtoFitNULL=TRUE",
  "--useSparseGRMforVarRatio=TRUE",
  "--isCateVarianceRatio=FALSE",
  "--LOCO=FALSE",
  "--IsOverwriteVarianceRatioFile=TRUE",
  paste0("--nThreads=", opt$n_threads)
)
message("Running: ", cmd)
status <- system(cmd)
if (status != 0) stop(sprintf("SAIGE step1 failed with status %s", status))

message(sprintf(
  "SAIGE null OK for %s: n=%d cases=%d controls=%d",
  opt$phenotype, nrow(pheno),
  sum(pheno[[opt$phenotype]] == 1), sum(pheno[[opt$phenotype]] == 0)
))
