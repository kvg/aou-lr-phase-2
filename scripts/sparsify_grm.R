#!/usr/bin/env Rscript
# Read PLINK2 square relatedness (.rel + .rel.id), sparsify, save RDS for GMMAT.
# Also emits kinship histogram / relationship-band counts / close-pair list for QC.

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list(
    rel = NA_character_,
    rel_id = NA_character_,
    samples = NA_character_,
    threshold = 0.05,
    close_pair_min = 0.2,
    out_rds = "grm_sparse.rds",
    out_counts = "grm_relatedness_summary.tsv",
    out_histogram = "grm_kinship_histogram.tsv",
    out_bands = "grm_relationship_bands.tsv",
    out_close_pairs = "grm_close_pairs.tsv"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    known <- c(
      "--rel", "--rel-id", "--samples", "--threshold", "--close-pair-min",
      "--out-rds", "--out-counts", "--out-histogram", "--out-bands",
      "--out-close-pairs"
    )
    if (key %in% known) {
      if (i == length(args)) stop(paste("Missing value for", key))
      val <- args[[i + 1]]
      if (key == "--rel") out$rel <- val
      if (key == "--rel-id") out$rel_id <- val
      if (key == "--samples") out$samples <- val
      if (key == "--threshold") out$threshold <- as.numeric(val)
      if (key == "--close-pair-min") out$close_pair_min <- as.numeric(val)
      if (key == "--out-rds") out$out_rds <- val
      if (key == "--out-counts") out$out_counts <- val
      if (key == "--out-histogram") out$out_histogram <- val
      if (key == "--out-bands") out$out_bands <- val
      if (key == "--out-close-pairs") out$out_close_pairs <- val
      i <- i + 2
    } else {
      stop(paste("Unknown arg:", key))
    }
  }
  if (is.na(out$rel) || is.na(out$rel_id) || is.na(out$samples)) {
    stop("Required: --rel --rel-id --samples")
  }
  out
}

# Hard-thresholding relatedness makes K indefinite. GMMAT glmmkin needs PD
# for Cholesky. A diagonal ridge K + λI is enough (λ > -λ_min), stays sparse,
# and is seconds at n~10k. Matrix::nearPD densifies and can take hours.
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

# Approximate relationship bands on GCTA / PLINK --make-rel scale (sibs ~0.5).
band_table <- function(kin) {
  breaks <- c(0.05, 0.088, 0.177, 0.354, 0.707, Inf)
  labels <- c(
    "cryptic_or_4th_plus_[0.05,0.088)",
    "approx_3rd_[0.088,0.177)",
    "approx_2nd_[0.177,0.354)",
    "approx_1st_[0.354,0.707)",
    "dup_or_MZ_[0.707,Inf)"
  )
  expected <- c(
    "distant / noise / 4th+",
    "~first cousins / 3rd",
    "~half-sib / avuncular / 2nd",
    "~full-sib / parent-offspring",
    "duplicate or MZ twin"
  )
  b <- cut(kin, breaks = breaks, right = FALSE, labels = labels)
  tab <- as.data.frame(table(band = b), stringsAsFactors = FALSE)
  names(tab)[2] <- "n_pairs"
  tab$expected_relationship <- expected[match(tab$band, labels)]
  tab$frac_of_pairs_ge_0.05 <- tab$n_pairs / max(1, length(kin))
  tab
}

opt <- parse_args(args)

ids_raw <- read.table(opt$rel_id, stringsAsFactors = FALSE, header = FALSE)
# PLINK2 --make-rel square writes either 1 or 2 columns (FID IID)
if (ncol(ids_raw) == 1) {
  rel_ids <- ids_raw[[1]]
} else {
  rel_ids <- ids_raw[[2]]
}
rel_ids <- as.character(rel_ids)

analysis <- scan(opt$samples, what = character(), quiet = TRUE)
analysis <- as.character(analysis)

missing <- setdiff(analysis, rel_ids)
if (length(missing) > 0) {
  stop(sprintf("%d analysis samples missing from GRM IDs (e.g. %s)",
               length(missing), paste(head(missing, 5), collapse = ",")))
}

message("Reading relatedness matrix: ", opt$rel)
K <- as.matrix(read.table(opt$rel, header = FALSE))
if (nrow(K) != length(rel_ids) || ncol(K) != length(rel_ids)) {
  stop(sprintf("GRM dim %dx%d does not match %d IDs", nrow(K), ncol(K), length(rel_ids)))
}
rownames(K) <- rel_ids
colnames(K) <- rel_ids

# Subset + reorder to analysis sample list
K <- K[analysis, analysis, drop = FALSE]

# Summary before sparsify
off <- K
diag(off) <- NA
n_pairs <- sum(!is.na(off)) / 2  # unique unordered pairs (upper+lower both filled)
# PLINK square rel is symmetric; count unique upper-tri pairs
ut <- upper.tri(off, diag = FALSE)
n_pairs_ut <- sum(ut)
vals_ge <- off[ut & !is.na(off) & off >= opt$threshold]
n_ge <- length(vals_ge)
max_off <- if (n_ge > 0) max(vals_ge) else max(off[ut], na.rm = TRUE)
diag_vals <- diag(K)
message(sprintf(
  "Off-diagonal pairs >= %.3f: %d / %d (max=%.4f)",
  opt$threshold, n_ge, n_pairs_ut, max_off
))

# Relationship-band counts among pairs >= threshold
bands <- band_table(vals_ge)
write.table(bands, opt$out_bands, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote bands: ", opt$out_bands)

# Histogram of kinship for pairs >= threshold
kin_max <- if (n_ge > 0) max(vals_ge) else opt$threshold
br <- seq(opt$threshold, kin_max, by = 0.025)
if (length(br) < 2 || tail(br, 1) <= kin_max) {
  br <- c(br, kin_max + 0.025)
}
br <- sort(unique(br))
h <- hist(vals_ge, breaks = br, plot = FALSE, right = FALSE, include.lowest = TRUE)
hist_df <- data.frame(
  kinship_bin_lo = h$breaks[-length(h$breaks)],
  kinship_bin_hi = h$breaks[-1],
  n_pairs = as.integer(h$counts)
)
write.table(hist_df, opt$out_histogram, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote histogram: ", opt$out_histogram)

# Close pairs for pedigree QC (default >= 0.2 ≈ strong 2nd / 1st)
close_idx <- which(ut & !is.na(off) & off >= opt$close_pair_min, arr.ind = TRUE)
if (nrow(close_idx) > 0) {
  close_df <- data.frame(
    id1 = rownames(off)[close_idx[, 1]],
    id2 = colnames(off)[close_idx[, 2]],
    kinship = off[close_idx],
    stringsAsFactors = FALSE
  )
  close_df <- close_df[order(-close_df$kinship), , drop = FALSE]
} else {
  close_df <- data.frame(
    id1 = character(),
    id2 = character(),
    kinship = numeric(),
    stringsAsFactors = FALSE
  )
}
write.table(close_df, opt$out_close_pairs, sep = "\t", quote = FALSE, row.names = FALSE)
message(sprintf(
  "Wrote close pairs (kinship >= %.3f): %s (%d pairs)",
  opt$close_pair_min, opt$out_close_pairs, nrow(close_df)
))

K[K < opt$threshold] <- 0
# Ensure diagonal is 1 for a kinship/relatedness matrix used as GMMAT kins.
# Thresholding can make the matrix non-PD; add the smallest diagonal ridge
# that restores a sparse Cholesky (no nearPD / densify).
diag(K) <- 1
K_sp <- Matrix::Matrix(K, sparse = TRUE)
K_sp <- regularize_kinship_for_gmmat(K_sp)
pd_ridge <- attr(K_sp, "pd_ridge")
if (is.null(pd_ridge)) pd_ridge <- NA_real_

saveRDS(K_sp, file = opt$out_rds)
message("Wrote sparse GRM: ", opt$out_rds)

summary_df <- data.frame(
  n_samples = length(analysis),
  threshold = opt$threshold,
  pd_ridge = pd_ridge,
  n_offdiag_pairs = n_pairs_ut,
  n_pairs_ge_threshold = n_ge,
  frac_pairs_ge_threshold = n_ge / max(1, n_pairs_ut),
  max_offdiag = max_off,
  mean_offdiag_ge_threshold = if (n_ge > 0) mean(vals_ge) else NA_real_,
  median_offdiag_ge_threshold = if (n_ge > 0) stats::median(vals_ge) else NA_real_,
  mean_diag = mean(diag_vals),
  min_diag = min(diag_vals),
  max_diag = max(diag_vals),
  n_close_pairs_ge_0.2 = nrow(close_df),
  stringsAsFactors = FALSE
)
write.table(summary_df, opt$out_counts, sep = "\t", quote = FALSE, row.names = FALSE)
message("Wrote summary: ", opt$out_counts)
