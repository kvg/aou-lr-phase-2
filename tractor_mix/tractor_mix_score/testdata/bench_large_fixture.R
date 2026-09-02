#!/usr/bin/env Rscript
# Larger synthetic fixture for R vs Rust timing comparison.

suppressPackageStartupMessages({
  library(GMMAT)
  library(Matrix)
  library(jsonlite)
})

out_dir <- if (length(commandArgs(trailingOnly = TRUE)) >= 1) {
  commandArgs(trailingOnly = TRUE)[[1]]
} else {
  file.path(tempdir(), "tractor_mix_bench_large")
}
args <- commandArgs(trailingOnly = TRUE)
n <- if (length(args) >= 2) as.integer(args[[2]]) else 800L
n_sites <- if (length(args) >= 3) as.integer(args[[3]]) else 2000L
n_anc <- if (length(args) >= 4) as.integer(args[[4]]) else 5L
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
    source = "gmmat", variance_ratio = 1,
    tractor_mix_score_sha = "bench-large"
  )
  write(toJSON(meta, auto_unbox = TRUE, pretty = TRUE), file.path(export_dir, "meta.json"))
  con <- file(file.path(export_dir, "sigma_i.csc.bin"), "wb")
  write_i32(con, n)
  write_i32(con, Matrix::nnzero(Sigma_i))
  write_i32(con, Sigma_i@p)
  write_i32(con, Sigma_i@i)
  write_f64(con, Sigma_i@x)
  close(con)
  con <- file(file.path(export_dir, "sigma_i_x.bin"), "wb")
  write_i32(con, n)
  write_i32(con, p)
  write_f64(con, Sigma_iX)
  close(con)
  con <- file(file.path(export_dir, "cov.bin"), "wb")
  write_i32(con, p)
  write_f64(con, cov_mat)
  close(con)
  con <- file(file.path(export_dir, "residuals.bin"), "wb")
  write_i32(con, n)
  write_f64(con, residuals)
  close(con)
}

set.seed(99)
samples <- sprintf("S%03d", seq_len(n))

K <- Matrix(0, n, n, sparse = TRUE)
for (i in seq_len(n)) {
  for (j in i:min(i + 8, n)) {
    if (i == j) {
      K[i, j] <- 1
    } else {
      v <- runif(1, 0, 0.08)
      K[i, j] <- v
      K[j, i] <- v
    }
  }
}
K <- K + Diagonal(n, 0.01)
dimnames(K) <- list(samples, samples)

sex <- rbinom(n, 1, 0.5)
age <- rnorm(n)
pheno <- data.frame(
  ID = samples,
  y = rbinom(n, 1, plogis(-0.5 + 0.4 * sex + 0.01 * age)),
  sex = sex,
  age = age
)

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

for (anc in seq_len(n_anc)) {
  con <- gzfile(file.path(dosage_dir, sprintf("anc_%02d.dosage.txt.gz", anc - 1)), "wt")
  on.exit(close(con), add = TRUE)
  writeLines(paste(header, collapse = "\t"), con)
  for (s in seq_len(n_sites)) {
    g <- sample(0:2, n, replace = TRUE)
    writeLines(
      paste(c("22", s * 100, paste0("v", s), "A", "G", as.character(g)), collapse = "\t"),
      con
    )
  }
  close(con)
  on.exit(NULL)
}

cat(
  "Fixture:", n, "samples,", n_sites, "sites x", n_anc, "ancestries;",
  "Sigma_i nnz:", Matrix::nnzero(Model_Null$Sigma_i), "\n"
)
