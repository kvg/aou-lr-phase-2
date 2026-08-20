#!/usr/bin/env Rscript
# Run vendored TractorMix.score on a fixture (oracle for parity tests).

args <- commandArgs(trailingOnly = TRUE)

`%||%` <- function(a, b) if (is.na(a)) b else a

parse_val <- function(flag) {
  i <- match(flag, args)
  if (is.na(i) || i == length(args)) NA_character_ else args[[i + 1]]
}

null_rds <- parse_val("--null-rds")
out_tsv <- parse_val("--out")
ac_threshold <- as.integer(parse_val("--ac-threshold") %||% "10")
n_core <- as.integer(parse_val("--n-core") %||% "1")
chunk_size <- as.integer(parse_val("--chunk-size") %||% "16")
oracle_r <- parse_val("--oracle-r")
if (is.na(oracle_r)) {
  oracle_r <- file.path(dirname(normalizePath(
    sub("--file=", "", grep("^--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1])
  )), "..", "oracle", "TractorMix.score.R")
}

dosage_files <- character(0)
idx <- which(args == "--dosage-files")
if (length(idx) == 1) {
  j <- idx + 1
  while (j <= length(args) && !startsWith(args[[j]], "--")) {
    dosage_files <- c(dosage_files, args[[j]])
    j <- j + 1
  }
}

if (is.na(null_rds) || length(dosage_files) < 1 || is.na(out_tsv)) {
  stop("Usage: run_oracle.R --null-rds X --dosage-files f1 f2 ... --out Y")
}

obj <- readRDS(null_rds)
source(oracle_r)
TractorMix.score(
  obj = obj,
  infiles = dosage_files,
  outfiles = out_tsv,
  AC_threshold = ac_threshold,
  n_core = n_core,
  chunk_size = chunk_size
)
