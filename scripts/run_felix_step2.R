#!/usr/bin/env Rscript
# Run FELIX step2_SPAtests.R on a FELIXla packed prefix (admixed / LAI tests).
#
# Keeps the native FELIX column set (P_cct_admixed_c, P_het_admixed_c,
# P_hom_admixed_c, per-ancestry *_c_*). Chromosome string must match the
# packed VCF contig (e.g. chr22, not 22).
#
# FELIX step2 only parallelizes when --idstoIncludeFile is set; a full-chrom
# FELIXla scan must use --nThreads=1 (the nThreads>1 branch otherwise no-ops).

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    felixla_prefix = NA_character_,
    chrom = "chr22",
    null_prefix = NA_character_,
    sample_file = NA_character_,
    sparse_grm = NA_character_,
    sparse_grm_ids = NA_character_,
    min_mac = 50L,
    n_ancestries = 5L,
    pvalcutoff_of_haplotype = 0.05,
    n_threads = 1L,
    step2_r = "/opt/SAIGE/extdata/step2_SPAtests.R",
    out_tsv = NA_character_,
    out_raw = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (i == length(args)) stop(paste("Missing value for", key))
    val <- args[[i + 1]]
    if (key == "--felixla-prefix") out$felixla_prefix <- val
    else if (key == "--chrom") out$chrom <- val
    else if (key == "--null-prefix") out$null_prefix <- val
    else if (key == "--sample-file") out$sample_file <- val
    else if (key == "--sparse-grm") out$sparse_grm <- val
    else if (key == "--sparse-grm-ids") out$sparse_grm_ids <- val
    else if (key == "--min-mac") out$min_mac <- as.integer(val)
    else if (key == "--n-ancestries") out$n_ancestries <- as.integer(val)
    else if (key == "--pvalcutoff-of-haplotype") {
      out$pvalcutoff_of_haplotype <- as.numeric(val)
    }
    else if (key == "--n-threads") out$n_threads <- as.integer(val)
    else if (key == "--step2-r") out$step2_r <- val
    else if (key == "--out-tsv") out$out_tsv <- val
    else if (key == "--out-raw") out$out_raw <- val
    else stop(paste("Unknown arg:", key))
    i <- i + 2
  }
  req <- c("felixla_prefix", "chrom", "null_prefix", "out_tsv",
           "sparse_grm", "sparse_grm_ids")
  for (r in req) {
    if (is.na(out[[r]]) || !nzchar(out[[r]])) stop(paste("Missing required", r))
  }
  if (is.na(out$out_raw) || !nzchar(out$out_raw)) {
    out$out_raw <- paste0(out$out_tsv, ".raw.txt")
  }
  out
}

opt <- parse_args(args)

if (!file.exists(opt$step2_r)) {
  fallback <- "/usr/local/bin/step2_SPAtests.R"
  if (file.exists(fallback)) {
    opt$step2_r <- fallback
  } else {
    stop(paste("step2_SPAtests.R not found:", opt$step2_r))
  }
}

meta <- paste0(opt$felixla_prefix, ".meta")
if (!file.exists(meta)) {
  stop(sprintf(
    "FELIXla prefix is missing %s (pass the prefix, not a filename)",
    meta
  ))
}

gmmat <- paste0(opt$null_prefix, ".rda")
if (!file.exists(gmmat)) {
  cands <- Sys.glob(paste0(opt$null_prefix, "*"))
  rda <- cands[grepl("\\.rda$", cands)]
  if (length(rda) < 1) {
    stop(paste("Null model .rda not found for prefix", opt$null_prefix))
  }
  gmmat <- rda[[1]]
}
var_ratio <- paste0(opt$null_prefix, ".varianceRatio.txt")
if (!file.exists(var_ratio)) {
  cands <- Sys.glob(paste0(opt$null_prefix, "*.varianceRatio.txt"))
  if (length(cands) < 1) stop("varianceRatio file not found")
  var_ratio <- cands[[1]]
}

if (!file.exists(opt$sparse_grm)) stop(paste("sparse GRM not found:", opt$sparse_grm))
if (!file.exists(opt$sparse_grm_ids)) {
  stop(paste("sparse GRM sample IDs not found:", opt$sparse_grm_ids))
}

# Full-chrom FELIXla scan: step2's nThreads>1 path requires --idstoIncludeFile
# and otherwise does no work. Pin the wrapper to 1 regardless of WDL cpu.
if (opt$n_threads != 1L) {
  message(sprintf(
    "Note: forcing FELIX step2 --nThreads=1 (requested %d); full-chrom FELIXla is single-threaded",
    opt$n_threads
  ))
}

subsample_arg <- ""
if (!is.na(opt$sample_file) && nzchar(opt$sample_file) && file.exists(opt$sample_file)) {
  subsample_arg <- paste0("--subSampleFile=", shQuote(opt$sample_file))
}

Sys.setenv(
  OMP_NUM_THREADS = "1",
  OPENBLAS_NUM_THREADS = "1",
  MKL_NUM_THREADS = "1"
)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(1)
  RhpcBLASctl::omp_set_num_threads(1)
}

dir.create(dirname(opt$out_raw), recursive = TRUE, showWarnings = FALSE)

cmd <- paste(
  "Rscript", shQuote(opt$step2_r),
  paste0("--FELIXlaPrefix=", shQuote(opt$felixla_prefix)),
  paste0("--chrom=", opt$chrom),
  "--is_admixed=TRUE",
  paste0("--number_of_ancestry=", opt$n_ancestries),
  paste0("--pvalcutoff_of_haplotype=", opt$pvalcutoff_of_haplotype),
  paste0("--minMAC=", opt$min_mac),
  paste0("--GMMATmodelFile=", shQuote(gmmat)),
  paste0("--varianceRatioFile=", shQuote(var_ratio)),
  paste0("--SAIGEOutputFile=", shQuote(opt$out_raw)),
  paste0("--sparseGRMFile=", shQuote(opt$sparse_grm)),
  paste0("--sparseGRMSampleIDFile=", shQuote(opt$sparse_grm_ids)),
  subsample_arg,
  "--nThreads=1",
  "--is_Firth_beta=TRUE",
  "--is_output_moreDetails=TRUE",
  "--LOCO=FALSE"
)
message("Running: ", cmd)
status <- system(cmd)
if (status != 0) stop(sprintf("FELIX step2 failed with status %s", status))

if (!file.exists(opt$out_raw)) {
  stop(paste("FELIX step2 did not write", opt$out_raw))
}

# Keep native FELIX columns for summarize (P_cct_admixed_c and *_c_*).
if (!file.exists(opt$out_tsv) || normalizePath(opt$out_raw) != normalizePath(opt$out_tsv, mustWork = FALSE)) {
  file.copy(opt$out_raw, opt$out_tsv, overwrite = TRUE)
}

hdr <- strsplit(readLines(opt$out_tsv, n = 1L), "\t", fixed = TRUE)[[1]]
need <- c("P_cct_admixed_c", "P_het_admixed_c", "P_hom_admixed_c")
missing <- setdiff(need, hdr)
if (length(missing) > 0) {
  stop(sprintf(
    "FELIX step2 output missing expected columns: %s (have: %s)",
    paste(missing, collapse = ","),
    paste(hdr, collapse = ",")
  ))
}

n <- as.integer(system(paste("wc -l <", shQuote(opt$out_tsv)), intern = TRUE))
message(sprintf(
  "Wrote FELIX results: %s (%d lines incl. header)",
  opt$out_tsv, n
))
