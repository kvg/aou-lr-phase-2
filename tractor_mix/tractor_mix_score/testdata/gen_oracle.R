#!/usr/bin/env Rscript
# Generate synthetic null + dosage fixtures for tractor-mix-score parity tests.

suppressPackageStartupMessages({
  library(GMMAT)
  library(Matrix)
  library(jsonlite)
})

args <- commandArgs(trailingOnly = TRUE)
out_dir <- if (length(args) >= 1) args[[1]] else {
  file.path(Sys.getenv("CARGO_MANIFEST_DIR", "."), "testdata", "fixture")
}
dir.create(out_dir, recursive = TRUE, show = FALSE)

write_f64 <- function(con, x) writeBin(as.double(x), con, size = 8, endian = "little")
write_i32 <- function(con, x) writeBin(as.integer(x), con, size = 4, endian = "little")

export_null_for_rust <- function(obj, export_dir) {
  dir.create(export_dir, recursive = TRUE, show = FALSE)
  n <- length(obj$id_include)
  p <- ncol(obj$X)
  Sigma_i <- as(obj$Sigma_i, "dgCMatrix")
  Sigma_iX <- as.matrix(obj$Sigma_iX)
  cov_mat <- as.matrix(obj$cov)
  residuals <- as.numeric(obj$scaled.residuals)
  writeLines(as.character(obj$id_include), file.path(export_dir, "id_include.txt"))
  meta <- list(
    n = n, p = p, nnz = Matrix::nnzero(Sigma_i), family = "binomial",
    tractor_mix_score_sha = "4adb8f1814d9315ecd7868eb729d52ec0c723719"
  )
  write(toJSON(meta, auto_unbox = TRUE, pretty = TRUE), file.path(export_dir, "meta.json"))
  con <- file(file.path(export_dir, "sigma_i.csc.bin"), "wb")
  write_i32(con, n); write_i32(con, Matrix::nnzero(Sigma_i))
  write_i32(con, Sigma_i@p); write_i32(con, Sigma_i@i); write_f64(con, Sigma_i@x)
  close(con)
  con <- file(file.path(export_dir, "sigma_i_x.bin"), "wb")
  write_i32(con, n); write_i32(con, p); write_f64(con, Sigma_iX)
  close(con)
  con <- file(file.path(export_dir, "cov.bin"), "wb")
  write_i32(con, p); write_f64(con, cov_mat)
  close(con)
  con <- file(file.path(export_dir, "residuals.bin"), "wb")
  write_i32(con, n); write_f64(con, residuals)
  close(con)
}

set.seed(42)
n <- 80
n_anc <- 5
n_sites <- 48
samples <- sprintf("S%03d", seq_len(n))

K <- Matrix(0, n, n, sparse = TRUE)
for (i in seq_len(n)) {
  for (j in i:min(i + 8, n)) {
    if (i == j) K[i, j] <- 1 else { v <- runif(1, 0, 0.08); K[i, j] <- v; K[j, i] <- v }
  }
}
K <- K + Diagonal(n, 0.01)
dimnames(K) <- list(samples, samples)

sex <- rbinom(n, 1, 0.5)
age <- rnorm(n)
eta <- -0.5 + 0.4 * sex + 0.01 * age
y <- rbinom(n, 1, plogis(eta))
pheno <- data.frame(ID = samples, y = y, sex = sex, age = age)

Model_Null <- glmmkin(
  fixed = y ~ sex + age,
  data = pheno,
  id = "ID",
  kins = K,
  family = binomial(link = "logit")
)
saveRDS(Model_Null, file.path(out_dir, "null.rds"))
export_null_for_rust(Model_Null, file.path(out_dir, "null_export"))

dosage_dir <- file.path(out_dir, "dosages")
dir.create(dosage_dir, show = FALSE)
header <- c("CHROM", "POS", "ID", "REF", "ALT", samples)

write_dosage <- function(path, rows) {
  con <- gzfile(path, "wt")
  on.exit(close(con))
  writeLines(paste(header, collapse = "\t"), con)
  for (i in seq_len(nrow(rows))) {
    writeLines(paste(c(rows$chrom[i], rows$pos[i], rows$id[i], rows$ref[i], rows$alt[i],
                       as.character(rows[i, samples])), collapse = "\t"), con)
  }
}

base <- data.frame(
  chrom = rep("22", n_sites),
  pos = as.character(seq(1000, by = 1000, length.out = n_sites)),
  id = sprintf("v%03d", seq_len(n_sites)),
  ref = "A",
  alt = "G",
  stringsAsFactors = FALSE
)

patterns <- list(
  rep(2L, n),
  c(rep(2L, 40), rep(0L, 40)),
  rep(0L, n),
  rep(1L, n)
)
for (i in seq_along(patterns)) base$id[i] <- sprintf("pattern%d", i)

for (anc in seq_len(n_anc)) {
  mat <- cbind(base, matrix(0L, n_sites, n))
  colnames(mat)[6:ncol(mat)] <- samples
  for (si in seq_len(n_sites)) {
    pat <- if (si <= length(patterns)) patterns[[si]] else sample(0:2, n, replace = TRUE)
    mat[si, samples] <- pat
  }
  write_dosage(
    file.path(dosage_dir, sprintf("anc_%02d.dosage.txt.gz", anc - 1)),
    mat
  )
}

cat("Wrote fixture to", out_dir, "\n")
