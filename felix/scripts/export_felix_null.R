#!/usr/bin/env Rscript
# Export a FELIX/SAIGE Step 1 modglmm (.rda) into the tractor-mix-score CSC layout.
#
# Writes:
#   id_include.txt, sigma_i.csc.bin, sigma_i_x.bin, cov.bin, residuals.bin,
#   xv.bin, xxvx_inv.bin, meta.json, variance_ratio.txt
#
# Sigma matches FELIX setSparseSigma_new (R/SAIGE_SPATest_Region.R):
#   binary:        Sigma = tau[2] * K + diag(1 / (mu*(1-mu)))
#   quantitative:  Sigma = tau[2] * K + diag(tau[1])
# Residuals are Step 1 residuals / tau[0] so Rust Score = gtilde' * residuals
# matches SAIGE S = gtilde' * res / tau[0]. One scalar variance ratio is stored
# (not MAC bins). XV / XXVX_inv come from obj.noK for getadjG.

suppressPackageStartupMessages({
  library(Matrix)
  library(data.table)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    null_rda = NA_character_,
    out_dir = "null_export",
    sparse_grm = NA_character_,
    sparse_grm_ids = NA_character_,
    variance_ratio = NA_character_,
    relatedness_cutoff = 0
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (i == length(args)) stop(paste("Missing value for", key))
    val <- args[[i + 1]]
    if (key == "--null-rda") out$null_rda <- val
    else if (key == "--out-dir") out$out_dir <- val
    else if (key == "--sparse-grm") out$sparse_grm <- val
    else if (key == "--sparse-grm-ids") out$sparse_grm_ids <- val
    else if (key == "--variance-ratio") out$variance_ratio <- val
    else if (key == "--relatedness-cutoff") out$relatedness_cutoff <- as.numeric(val)
    else if (key == "--maxiter-pcg" || key == "--tol-pcg") {
      # kept for CLI compatibility; inversion uses direct solve()
    } else {
      stop(paste("Unknown arg:", key))
    }
    i <- i + 2
  }
  if (is.na(out$null_rda) || !nzchar(out$null_rda)) {
    stop("Required: --null-rda")
  }
  out
}

write_f64 <- function(con, x) {
  writeBin(as.double(x), con, size = 8, endian = "little")
}

write_i32 <- function(con, x) {
  writeBin(as.integer(x), con, size = 4, endian = "little")
}

parse_scalar_variance_ratio <- function(path) {
  if (is.na(path) || !nzchar(path) || !file.exists(path)) {
    return(1)
  }
  dt <- tryCatch(
    data.table::fread(path, header = FALSE, data.table = FALSE),
    error = function(e) NULL
  )
  if (is.null(dt) || nrow(dt) < 1) {
    return(1)
  }
  if (ncol(dt) >= 2) {
    typ <- as.character(dt[[2]])
    sparse <- which(typ == "sparse")
    if (length(sparse) > 0) {
      return(as.numeric(dt[[1]][sparse[[1]]]))
    }
    nul <- which(typ == "null")
    if (length(nul) > 0) {
      return(as.numeric(dt[[1]][nul[[1]]]))
    }
  }
  as.numeric(dt[[1]][[1]])
}

read_sparse_grm <- function(mtx_path, ids_path, sample_ids, relatedness_cutoff = 0) {
  if (!file.exists(mtx_path)) stop("sparse GRM not found: ", mtx_path)
  if (!file.exists(ids_path)) stop("sparse GRM IDs not found: ", ids_path)
  grm_ids <- scan(ids_path, what = character(), quiet = TRUE)
  sparse_grm <- Matrix::readMM(mtx_path)
  sparse_grm <- as(sparse_grm, "dgCMatrix")
  if (length(grm_ids) != nrow(sparse_grm)) {
    stop(sprintf("GRM ID count %d != mtx n=%d", length(grm_ids), nrow(sparse_grm)))
  }
  idx <- match(sample_ids, grm_ids)
  if (any(is.na(idx))) {
    stop(sprintf("%d samples missing from sparse GRM", sum(is.na(idx))))
  }
  sparse_grm <- as(sparse_grm[idx, idx, drop = FALSE], "dgTMatrix")
  if (relatedness_cutoff > 0) {
    drop <- which(sparse_grm@x < relatedness_cutoff)
    if (length(drop) > 0) {
      message(sprintf(
        "Removing %d sparse-GRM elements < %g",
        length(drop), relatedness_cutoff
      ))
      sparse_grm@x <- sparse_grm@x[-drop]
      sparse_grm@i <- sparse_grm@i[-drop]
      sparse_grm@j <- sparse_grm@j[-drop]
    }
  }
  as(sparse_grm, "dgCMatrix")
}

reconstruct_nok <- function(mod, mu2) {
  X <- as.matrix(mod$X)
  v <- as.vector(mu2)
  xv <- t(X * v)
  xvx <- t(X) %*% t(xv)
  xxvx_inv <- X %*% solve(xvx)
  list(XV = xv, XXVX_inv = xxvx_inv)
}

build_export_obj <- function(mod, sparse_grm = NA_character_, sparse_grm_ids = NA_character_,
                             relatedness_cutoff = 0) {
  n <- length(mod$sampleID)
  ids <- as.character(mod$sampleID)
  X <- as.matrix(mod$X)
  p <- ncol(X)
  tau <- as.numeric(mod$theta)
  if (length(tau) < 1) stop("modglmm$theta is empty")
  tau0 <- tau[[1]]
  tau_k <- if (length(tau) >= 2) tau[[2]] else 0
  mu <- as.vector(mod$fitted.values)
  if (length(mu) != n) {
    stop(sprintf("fitted.values length %d != n %d", length(mu), n))
  }

  if (identical(mod$traitType, "binary")) {
    mu2 <- pmax(mu * (1 - mu), 1e-8)
    family <- "binomial"
    trait_type <- "binary"
    d_add <- 1 / mu2
  } else if (identical(mod$traitType, "quantitative")) {
    mu2 <- rep(1 / pmax(tau0, 1e-8), n)
    family <- "gaussian"
    trait_type <- "quantitative"
    d_add <- rep(tau0, n)
  } else {
    stop(sprintf("unsupported traitType %s (binary or quantitative)", mod$traitType))
  }

  use_sparse <- !is.na(sparse_grm) && nzchar(sparse_grm) &&
    !is.na(sparse_grm_ids) && nzchar(sparse_grm_ids)

  if (use_sparse) {
    message(sprintf("Building Sigma from sparse GRM (n=%d, tau_k=%g)", n, tau_k))
    K <- read_sparse_grm(sparse_grm, sparse_grm_ids, ids, relatedness_cutoff)
    Sigma <- tau_k * K
    diag(Sigma) <- diag(Sigma) + d_add
  } else {
    message("No sparse GRM provided; Sigma is diagonal working-weight / residual variance")
    Sigma <- Diagonal(n, x = d_add)
  }

  message(sprintf("Inverting %d x %d Sigma for CSC export", n, n))
  Sigma_i <- tryCatch(
    solve(Sigma),
    error = function(e) {
      message("Sigma solve failed; retrying with 1e-8 diagonal ridge")
      solve(Sigma + Diagonal(n, x = 1e-8))
    }
  )
  Sigma_i <- as(Sigma_i, "dgCMatrix")
  Sigma_iX <- as.matrix(Sigma_i %*% X)
  xtsx <- crossprod(X, Sigma_iX)
  cov_mat <- tryCatch(
    solve(xtsx),
    error = function(e) {
      message("X' Sigma_i X solve failed; using ginv-style ridge")
      solve(xtsx + diag(1e-8, p))
    }
  )

  residuals <- as.numeric(mod$residuals)
  if (length(residuals) != n) {
    stop(sprintf("residuals length %d != n %d", length(residuals), n))
  }
  residuals <- residuals / pmax(tau0, 1e-12)

  if (!is.null(mod$obj.noK$XV) && !is.null(mod$obj.noK$XXVX_inv)) {
    xv <- as.matrix(mod$obj.noK$XV)
    xxvx_inv <- as.matrix(mod$obj.noK$XXVX_inv)
  } else {
    message("obj.noK missing XV/XXVX_inv; reconstructing getadjG matrices")
    nok <- reconstruct_nok(mod, mu2)
    xv <- nok$XV
    xxvx_inv <- nok$XXVX_inv
  }
  if (nrow(xv) != p || ncol(xv) != n) {
    stop(sprintf("XV dims %dx%d != p=%d n=%d", nrow(xv), ncol(xv), p, n))
  }
  if (nrow(xxvx_inv) != n || ncol(xxvx_inv) != p) {
    stop(sprintf("XXVX_inv dims %dx%d != n=%d p=%d", nrow(xxvx_inv), ncol(xxvx_inv), n, p))
  }

  list(
    family = family,
    trait_type = trait_type,
    id_include = ids,
    Sigma_i = Sigma_i,
    Sigma_iX = Sigma_iX,
    cov = as.matrix(cov_mat),
    residuals = residuals,
    xv = xv,
    xxvx_inv = xxvx_inv,
    n = n,
    p = p,
    tau0 = tau0
  )
}

export_null_for_rust <- function(obj, out_dir, variance_ratio = 1) {
  if (!is.finite(variance_ratio) || variance_ratio <= 0) {
    stop("variance_ratio must be a positive finite scalar")
  }
  dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

  n <- obj$n
  p <- obj$p
  Sigma_i <- obj$Sigma_i

  writeLines(obj$id_include, file.path(out_dir, "id_include.txt"))

  meta <- list(
    n = n,
    p = p,
    nnz = Matrix::nnzero(Sigma_i),
    family = obj$family,
    trait_type = obj$trait_type,
    source = "felix",
    variance_ratio = variance_ratio,
    tau0 = obj$tau0,
    tractor_mix_score_sha = "felix-export-v1"
  )
  write(toJSON(meta, auto_unbox = TRUE, pretty = TRUE), file.path(out_dir, "meta.json"))

  con <- file(file.path(out_dir, "sigma_i.csc.bin"), "wb")
  write_i32(con, n)
  write_i32(con, Matrix::nnzero(Sigma_i))
  write_i32(con, Sigma_i@p)
  write_i32(con, Sigma_i@i)
  write_f64(con, Sigma_i@x)
  close(con)

  con <- file(file.path(out_dir, "sigma_i_x.bin"), "wb")
  write_i32(con, n)
  write_i32(con, p)
  write_f64(con, obj$Sigma_iX)
  close(con)

  con <- file(file.path(out_dir, "cov.bin"), "wb")
  write_i32(con, p)
  write_f64(con, obj$cov)
  close(con)

  con <- file(file.path(out_dir, "residuals.bin"), "wb")
  write_i32(con, n)
  write_f64(con, obj$residuals)
  close(con)

  con <- file(file.path(out_dir, "xv.bin"), "wb")
  write_i32(con, p)
  write_i32(con, n)
  write_f64(con, obj$xv)
  close(con)

  con <- file(file.path(out_dir, "xxvx_inv.bin"), "wb")
  write_i32(con, n)
  write_i32(con, p)
  write_f64(con, obj$xxvx_inv)
  close(con)

  writeLines(format(variance_ratio, scientific = TRUE, digits = 16),
             file.path(out_dir, "variance_ratio.txt"))

  message(sprintf(
    "Wrote FELIX null export: %s (n=%d, p=%d, nnz=%d, family=%s, var_ratio=%g)",
    out_dir, n, p, meta$nnz, obj$family, variance_ratio
  ))
}

opt <- parse_args(args)
load(opt$null_rda)
if (!exists("modglmm")) stop("Expected object modglmm in ", opt$null_rda)

mod <- build_export_obj(
  modglmm,
  sparse_grm = opt$sparse_grm,
  sparse_grm_ids = opt$sparse_grm_ids,
  relatedness_cutoff = opt$relatedness_cutoff
)
vr_path <- opt$variance_ratio
if (is.na(vr_path) || !nzchar(vr_path)) {
  cand <- sub("\\.rda$", ".varianceRatio.txt", opt$null_rda)
  if (file.exists(cand)) vr_path <- cand
}
vr <- parse_scalar_variance_ratio(vr_path)
export_null_for_rust(mod, opt$out_dir, variance_ratio = vr)
