#!/usr/bin/env Rscript
# Fit GMMAT null model (binomial) and export sparse Sigma_i artifacts for tractor-mix-score.
# FELIX/SAIGE Step 1 .rda: pass --step1-rda (delegates to export_felix_null.R).

suppressPackageStartupMessages({
  library(GMMAT)
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    pheno_cov = NA_character_,
    phenotype = NA_character_,
    covariates = NA_character_,
    grm_rds = NA_character_,
    out_null_rds = "null_model.rds",
    out_null_export = "null_export",
    step1_rda = NA_character_,
    sparse_grm = NA_character_,
    sparse_grm_ids = NA_character_,
    variance_ratio = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (i == length(args)) stop(paste("Missing value for", key))
    val <- args[[i + 1]]
    if (key == "--pheno-cov") out$pheno_cov <- val
    else if (key == "--phenotype") out$phenotype <- val
    else if (key == "--covariates") out$covariates <- val
    else if (key == "--grm-rds") out$grm_rds <- val
    else if (key == "--out-null-rds") out$out_null_rds <- val
    else if (key == "--out-null-export") out$out_null_export <- val
    else if (key == "--step1-rda") out$step1_rda <- val
    else if (key == "--sparse-grm") out$sparse_grm <- val
    else if (key == "--sparse-grm-ids") out$sparse_grm_ids <- val
    else if (key == "--variance-ratio") out$variance_ratio <- val
    else stop(paste("Unknown arg:", key))
    i <- i + 2
  }
  if (!is.na(out$step1_rda) && nzchar(out$step1_rda)) {
    return(out)
  }
  if (is.na(out$pheno_cov) || is.na(out$phenotype) || is.na(out$covariates) ||
      is.na(out$grm_rds)) {
    stop("Required: --pheno-cov --phenotype --covariates --grm-rds (or --step1-rda)")
  }
  out
}

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

write_f64 <- function(con, x) {
  writeBin(as.double(x), con, size = 8, endian = "little")
}

write_i32 <- function(con, x) {
  writeBin(as.integer(x), con, size = 4, endian = "little")
}

export_null_for_rust <- function(obj, out_dir) {
  if (!is.null(obj$P)) {
    stop("Dense GRM null (obj$P non-NULL) is not supported by tractor-mix-score")
  }
  if (!any(grepl("binomial", as.character(obj$call)))) {
    stop("Only binomial null models are supported by tractor-mix-score")
  }

  dir.create(out_dir, recursive = TRUE, show = FALSE)

  n <- length(obj$id_include)
  p <- ncol(obj$X)
  Sigma_i <- as(obj$Sigma_i, "dgCMatrix")
  Sigma_iX <- as.matrix(obj$Sigma_iX)
  cov_mat <- as.matrix(obj$cov)
  residuals <- as.numeric(obj$scaled.residuals)

  if (length(residuals) != n) {
    stop(sprintf("scaled.residuals length %d != n %d", length(residuals), n))
  }
  if (nrow(Sigma_iX) != n || ncol(Sigma_iX) != p) {
    stop(sprintf("Sigma_iX dims %dx%d != n=%d p=%d", nrow(Sigma_iX), ncol(Sigma_iX), n, p))
  }
  if (nrow(cov_mat) != p || ncol(cov_mat) != p) {
    stop(sprintf("cov dims %dx%d != p=%d", nrow(cov_mat), ncol(cov_mat), p))
  }

  writeLines(as.character(obj$id_include), file.path(out_dir, "id_include.txt"))

  meta <- list(
    n = n,
    p = p,
    nnz = Matrix::nnzero(Sigma_i),
    family = "binomial",
    source = "gmmat",
    variance_ratio = 1,
    tractor_mix_score_sha = "4adb8f1814d9315ecd7868eb729d52ec0c723719"
  )
  write(toJSON(meta, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "meta.json"))
  writeLines("1", file.path(out_dir, "variance_ratio.txt"))

  # CSC binary: n, nnz, colptr (0-based), rowidx (0-based), values
  con <- file(file.path(out_dir, "sigma_i.csc.bin"), "wb")
  on.exit(close(con), add = TRUE)
  write_i32(con, n)
  write_i32(con, Matrix::nnzero(Sigma_i))
  write_i32(con, Sigma_i@p)
  write_i32(con, Sigma_i@i)
  write_f64(con, Sigma_i@x)
  close(con)
  on.exit(NULL)

  con <- file(file.path(out_dir, "sigma_i_x.bin"), "wb")
  on.exit(close(con), add = TRUE)
  write_i32(con, n)
  write_i32(con, p)
  write_f64(con, Sigma_iX)
  close(con)
  on.exit(NULL)

  con <- file(file.path(out_dir, "cov.bin"), "wb")
  on.exit(close(con), add = TRUE)
  write_i32(con, p)
  write_f64(con, cov_mat)
  close(con)
  on.exit(NULL)

  con <- file(file.path(out_dir, "residuals.bin"), "wb")
  on.exit(close(con), add = TRUE)
  write_i32(con, n)
  write_f64(con, residuals)
  close(con)
  on.exit(NULL)

  message(sprintf("Wrote null export: %s (n=%d, p=%d, nnz=%d)", out_dir, n, p, meta$nnz))
}

opt <- parse_args(args)

if (!is.na(opt$step1_rda) && nzchar(opt$step1_rda)) {
  script_file <- sub("^--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  script_dir <- if (length(script_file) && nzchar(script_file) && !is.na(script_file)) {
    dirname(normalizePath(script_file, mustWork = FALSE))
  } else {
    getwd()
  }
  export_r <- file.path(script_dir, "export_felix_null.R")
  if (!file.exists(export_r)) {
    export_r <- file.path(dirname(script_dir), "felix", "scripts", "export_felix_null.R")
  }
  if (!file.exists(export_r)) {
    export_r <- file.path("/opt/felix_scripts/export_felix_null.R")
  }
  if (!file.exists(export_r)) {
    stop("export_felix_null.R not found next to fit_null.R")
  }
  cmd <- paste(
    "Rscript", shQuote(export_r),
    "--null-rda", shQuote(opt$step1_rda),
    "--out-dir", shQuote(opt$out_null_export)
  )
  if (!is.na(opt$sparse_grm) && nzchar(opt$sparse_grm)) {
    cmd <- paste(cmd, "--sparse-grm", shQuote(opt$sparse_grm))
  }
  if (!is.na(opt$sparse_grm_ids) && nzchar(opt$sparse_grm_ids)) {
    cmd <- paste(cmd, "--sparse-grm-ids", shQuote(opt$sparse_grm_ids))
  }
  if (!is.na(opt$variance_ratio) && nzchar(opt$variance_ratio)) {
    cmd <- paste(cmd, "--variance-ratio", shQuote(opt$variance_ratio))
  }
  message("Delegating FELIX Step 1 export: ", cmd)
  status <- system(cmd)
  if (status != 0) stop(sprintf("export_felix_null.R failed with status %s", status))
  quit(save = "no", status = 0)
}

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

export_null_for_rust(Model_Null, opt$out_null_export)
