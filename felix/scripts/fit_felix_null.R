#!/usr/bin/env Rscript
# Fit a FELIX / SAIGE null model (step 1) for one phenotype (binary or quantitative).
#
# Aligns pheno_cov to analysis_samples, runs step1_fitNULLGLMM.R, then exports
# null_export artifacts for tractor-mix-score via export_felix_null.R.

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
    trait_type = "binary",
    step1_r = "step1_fitNULLGLMM.R",
    export_r = NA_character_,
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
    else if (key == "--trait-type") out$trait_type <- val
    else if (key == "--step1-r") out$step1_r <- val
    else if (key == "--export-r") out$export_r <- val
    else if (key == "--n-threads") out$n_threads <- as.integer(val)
    else if (key == "--out-prefix") out$out_prefix <- val
    else stop(paste("Unknown arg:", key))
    i <- i + 2
  }
  req <- c(
    "pheno_cov", "phenotype", "covariates", "analysis_samples",
    "plink_prefix", "sparse_grm", "sparse_grm_ids", "out_prefix"
  )
  for (r in req) {
    if (is.na(out[[r]]) || !nzchar(out[[r]])) stop(paste("Missing required", r))
  }
  if (!(out$trait_type %in% c("binary", "quantitative"))) {
    stop("--trait-type must be binary or quantitative")
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

if (opt$trait_type == "binary") {
  vals <- sort(unique(pheno[[opt$phenotype]]))
  if (!all(vals %in% c(0, 1))) {
    stop(sprintf(
      "Phenotype %s is not binary 0/1 after filtering: %s",
      opt$phenotype, paste(vals, collapse = ",")
    ))
  }
}

dir.create(dirname(opt$out_prefix), recursive = TRUE, showWarnings = FALSE)
pheno_path <- paste0(opt$out_prefix, ".pheno.tsv")
fwrite(pheno[, keep_cols, drop = FALSE], pheno_path, sep = "\t", na = "NA", quote = FALSE)

sample_path <- paste0(opt$out_prefix, ".samples.txt")
writeLines(pheno$ID, sample_path)

covar_list <- paste(covars, collapse = ",")
meta_path <- paste0(opt$out_prefix, ".null_meta.tsv")
n_cases <- if (opt$trait_type == "binary") sum(pheno[[opt$phenotype]] == 1) else NA_integer_
n_controls <- if (opt$trait_type == "binary") sum(pheno[[opt$phenotype]] == 0) else NA_integer_
fwrite(
  data.frame(
    phenotype = opt$phenotype,
    trait_type = opt$trait_type,
    n_samples = nrow(pheno),
    n_cases = n_cases,
    n_controls = n_controls,
    covariates = covar_list,
    stringsAsFactors = FALSE
  ),
  meta_path,
  sep = "\t",
  quote = FALSE
)

saige_trait <- if (opt$trait_type == "binary") "binary" else "quantitative"
cmd <- paste(
  "Rscript", shQuote(opt$step1_r),
  paste0("--plinkFile=", shQuote(opt$plink_prefix)),
  paste0("--phenoFile=", shQuote(pheno_path)),
  paste0("--phenoCol=", opt$phenotype),
  paste0("--covarColList=", covar_list),
  "--sampleIDColinphenoFile=ID",
  paste0("--traitType=", saige_trait),
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
if (status != 0) stop(sprintf("FELIX/SAIGE step1 failed with status %s", status))

export_r <- opt$export_r
if (is.na(export_r) || !nzchar(export_r)) {
  script_dir <- dirname(normalizePath(
    sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]),
    mustWork = FALSE
  ))
  export_r <- file.path(script_dir, "export_felix_null.R")
  if (!file.exists(export_r)) {
    export_r <- file.path("/opt/felix_scripts/export_felix_null.R")
  }
}
null_export_dir <- paste0(opt$out_prefix, ".null_export")
export_cmd <- paste(
  "Rscript", shQuote(export_r),
  "--null-rda", shQuote(paste0(opt$out_prefix, ".rda")),
  "--out-dir", shQuote(null_export_dir),
  "--sparse-grm", shQuote(opt$sparse_grm),
  "--sparse-grm-ids", shQuote(opt$sparse_grm_ids),
  "--variance-ratio", shQuote(paste0(opt$out_prefix, ".varianceRatio.txt"))
)
message("Running: ", export_cmd)
status <- system(export_cmd)
if (status != 0) stop(sprintf("export_felix_null failed with status %s", status))

message(sprintf(
  "FELIX null OK for %s (%s): n=%d",
  opt$phenotype, opt$trait_type, nrow(pheno)
))
