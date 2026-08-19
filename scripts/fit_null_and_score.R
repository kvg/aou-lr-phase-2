#!/usr/bin/env Rscript
# Fit GMMAT null model and run unconditional TractorMix.score for one phenotype.

suppressPackageStartupMessages({
  library(GMMAT)
  library(Matrix)
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    pheno_cov = NA_character_,
    phenotype = NA_character_,
    covariates = NA_character_,
    grm_rds = NA_character_,
    dosage_files = character(0),
    tractor_mix_score_r = "/opt/Tractor-Mix/TractorMix.score.R",
    ac_threshold = 50L,
    n_core = 4L,
    out_tsv = "tractor_mix_results.tsv",
    out_null_rds = "null_model.rds"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (key == "--dosage-files") {
      # remaining until next -- or end
      i <- i + 1
      while (i <= length(args) && !startsWith(args[[i]], "--")) {
        out$dosage_files <- c(out$dosage_files, args[[i]])
        i <- i + 1
      }
      next
    }
    if (i == length(args)) stop(paste("Missing value for", key))
    val <- args[[i + 1]]
    if (key == "--pheno-cov") out$pheno_cov <- val
    else if (key == "--phenotype") out$phenotype <- val
    else if (key == "--covariates") out$covariates <- val
    else if (key == "--grm-rds") out$grm_rds <- val
    else if (key == "--tractor-mix-score-r") out$tractor_mix_score_r <- val
    else if (key == "--ac-threshold") out$ac_threshold <- as.integer(val)
    else if (key == "--n-core") out$n_core <- as.integer(val)
    else if (key == "--out-tsv") out$out_tsv <- val
    else if (key == "--out-null-rds") out$out_null_rds <- val
    else stop(paste("Unknown arg:", key))
    i <- i + 2
  }
  if (is.na(out$pheno_cov) || is.na(out$phenotype) || is.na(out$covariates) ||
      is.na(out$grm_rds) || length(out$dosage_files) < 1) {
    stop("Required: --pheno-cov --phenotype --covariates --grm-rds --dosage-files ...")
  }
  out
}

# Thresholded PLINK relatedness is often indefinite. GMMAT glmmkin needs a PD
# kinship matrix for Cholesky. A diagonal ridge K + λI with λ > -λ_min is
# enough, keeps the sparsity pattern, and is seconds at n~10k. Do not use
# Matrix::nearPD here: it densifies and can take hours.
kinship_is_pd <- function(K) {
  tryCatch({
    Matrix::chol(Matrix::forceSymmetric(as(K, "dsCMatrix")))
    TRUE
  }, error = function(e) FALSE)
}

find_pd_ridge <- function(K, max_ridge = 50) {
  n <- nrow(K)
  if (kinship_is_pd(K)) {
    return(0)
  }
  lo <- 0
  hi <- 1e-4
  while (!kinship_is_pd(K + Matrix::Diagonal(n, hi))) {
    lo <- hi
    hi <- hi * 2
    if (lo >= max_ridge) {
      stop(sprintf(
        "GRM not PD after diagonal ridge %.4g (n=%d). Check the relatedness matrix.",
        max_ridge, n
      ))
    }
  }
  hi <- min(hi, max_ridge)
  for (i in seq_len(16)) {
    mid <- (lo + hi) / 2
    if (kinship_is_pd(K + Matrix::Diagonal(n, mid))) {
      hi <- mid
    } else {
      lo <- mid
    }
  }
  hi
}

regularize_kinship_for_gmmat <- function(K) {
  K <- as(K, "dgCMatrix")
  K <- (K + Matrix::t(K)) / 2
  dn <- dimnames(K)
  n <- nrow(K)
  diag(K) <- pmax(Matrix::diag(K), 1)
  ridge <- find_pd_ridge(K)
  if (ridge > 0) {
    K <- K + Matrix::Diagonal(n, ridge)
  }
  nnz <- Matrix::nnzero(K)
  if (ridge == 0) {
    message(sprintf("GRM is positive definite (no ridge; n=%d, nnz=%d)", n, nnz))
  } else {
    message(sprintf(
      "Applied GRM diagonal ridge %.4g for positive definiteness (n=%d, nnz=%d); sparsity preserved",
      ridge, n, nnz
    ))
  }
  out <- as(K, "sparseMatrix")
  dimnames(out) <- dn
  attr(out, "pd_ridge") <- ridge
  out
}

opt <- parse_args(args)

covars <- scan(opt$covariates, what = character(), quiet = TRUE)
pheno <- fread(opt$pheno_cov, sep = "\t", header = TRUE, data.table = FALSE)
if (!("ID" %in% names(pheno))) stop("pheno_cov must contain ID column")
if (!(opt$phenotype %in% names(pheno))) {
  stop(sprintf("Phenotype %s not in pheno_cov columns", opt$phenotype))
}
missing_cov <- setdiff(covars, names(pheno))
if (length(missing_cov) > 0) {
  stop(sprintf("Missing covariate columns: %s", paste(missing_cov, collapse = ",")))
}

# Complete cases for this phenotype + covariates
keep_cols <- c("ID", opt$phenotype, covars)
cc <- stats::complete.cases(pheno[, keep_cols])
pheno <- pheno[cc, , drop = FALSE]
pheno$ID <- as.character(pheno$ID)
pheno[[opt$phenotype]] <- as.numeric(pheno[[opt$phenotype]])
for (cv in covars) {
  pheno[[cv]] <- as.numeric(pheno[[cv]])
}

message(sprintf(
  "Phenotype %s: n=%d cases=%d controls=%d",
  opt$phenotype, nrow(pheno),
  sum(pheno[[opt$phenotype]] == 1, na.rm = TRUE),
  sum(pheno[[opt$phenotype]] == 0, na.rm = TRUE)
))

GRM <- readRDS(opt$grm_rds)
grm_ids <- rownames(GRM)
if (is.null(grm_ids)) stop("GRM RDS must have rownames = sample IDs")

missing_grm <- setdiff(pheno$ID, grm_ids)
if (length(missing_grm) > 0) {
  stop(sprintf("%d phenotype samples missing from GRM", length(missing_grm)))
}

# Align GRM to phenotype sample order and ensure PD for glmmkin
GRM <- GRM[pheno$ID, pheno$ID, drop = FALSE]
GRM <- regularize_kinship_for_gmmat(GRM)

fixed_formula <- as.formula(
  paste(opt$phenotype, "~", paste(covars, collapse = " + "))
)
message("Null model: ", deparse(fixed_formula))

Model_Null <- glmmkin(
  fixed = fixed_formula,
  data = pheno,
  id = "ID",
  kins = GRM,
  family = binomial(link = "logit")
)
saveRDS(Model_Null, opt$out_null_rds)
message("Wrote null model: ", opt$out_null_rds)

source(opt$tractor_mix_score_r)
# Decompress .gz dosages to plain txt if needed (TractorMix fread handles gz, but
# ensure paths are passed as-is; TractorMix.score supports .gz via zcat).
TractorMix.score(
  obj = Model_Null,
  infiles = opt$dosage_files,
  outfiles = opt$out_tsv,
  AC_threshold = opt$ac_threshold,
  n_core = opt$n_core
)
message("Wrote Tractor-Mix results: ", opt$out_tsv)
