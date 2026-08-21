#!/usr/bin/env Rscript
# Run SAIGE step 2 on a VCF and normalize association columns for QC.

suppressPackageStartupMessages({
  library(data.table)
})

args <- commandArgs(trailingOnly = TRUE)

parse_args <- function(args) {
  out <- list(
    vcf = NA_character_,
    chrom = "chr22",
    null_prefix = NA_character_,
    sample_file = NA_character_,
    sparse_grm = NA_character_,
    sparse_grm_ids = NA_character_,
    min_mac = 20L,
    n_threads = 8L,
    step2_r = "/opt/SAIGE/extdata/step2_SPAtests.R",
    out_tsv = NA_character_,
    out_raw = NA_character_
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    if (i == length(args)) stop(paste("Missing value for", key))
    val <- args[[i + 1]]
    if (key == "--vcf") out$vcf <- val
    else if (key == "--chrom") out$chrom <- val
    else if (key == "--null-prefix") out$null_prefix <- val
    else if (key == "--sample-file") out$sample_file <- val
    else if (key == "--sparse-grm") out$sparse_grm <- val
    else if (key == "--sparse-grm-ids") out$sparse_grm_ids <- val
    else if (key == "--min-mac") out$min_mac <- as.integer(val)
    else if (key == "--n-threads") out$n_threads <- as.integer(val)
    else if (key == "--step2-r") out$step2_r <- val
    else if (key == "--out-tsv") out$out_tsv <- val
    else if (key == "--out-raw") out$out_raw <- val
    else stop(paste("Unknown arg:", key))
    i <- i + 2
  }
  req <- c("vcf", "null_prefix", "out_tsv", "sparse_grm", "sparse_grm_ids")
  for (r in req) {
    if (is.na(out[[r]]) || !nzchar(out[[r]])) stop(paste("Missing required", r))
  }
  if (is.na(out$out_raw) || !nzchar(out$out_raw)) {
    out$out_raw <- paste0(out$out_tsv, ".raw.txt")
  }
  out
}

pick_col <- function(df, candidates) {
  for (c in candidates) {
    if (c %in% names(df)) return(c)
  }
  NA_character_
}

opt <- parse_args(args)

# FLARE anc VCFs often declare ##FORMAT=<ID=GT,...> twice; SAIGE/htslib can fail.
vcf_in <- opt$vcf
dir.create("saige_vcf_fix", showWarnings = FALSE)
hdr <- "saige_vcf_fix/header.fixed.txt"
vcf_fixed <- "saige_vcf_fix/plink_in.vcf.gz"
status <- system(paste(
  "bcftools view -h", shQuote(vcf_in),
  "| awk '/^##FORMAT=<ID=GT,/{ if (gt++) next } { print }' >", shQuote(hdr)
))
if (status != 0) stop("bcftools view -h / awk header fix failed")
n_orig <- as.integer(system(paste(
  "bcftools view -h", shQuote(vcf_in), "| grep -c '^##FORMAT=<ID=GT,' || true"
), intern = TRUE)[[1]])
n_fix <- as.integer(system(paste(
  "grep -c '^##FORMAT=<ID=GT,'", shQuote(hdr), "|| true"
), intern = TRUE)[[1]])
message(sprintf("FORMAT/GT header lines: %d -> %d", n_orig, n_fix))
status <- system(paste("bcftools reheader -h", shQuote(hdr), "-o", shQuote(vcf_fixed), shQuote(vcf_in)))
if (status != 0) stop("bcftools reheader failed")
status <- system(paste("bcftools index -t", shQuote(vcf_fixed)))
if (status != 0) {
  status <- system(paste("bcftools index -c", shQuote(vcf_fixed)))
  if (status != 0) stop("bcftools index failed on fixed VCF")
}
opt$vcf <- vcf_fixed

# Ensure VCF index exists beside the localized file
idx_tbi <- paste0(opt$vcf, ".tbi")
idx_csi <- paste0(opt$vcf, ".csi")
if (!file.exists(idx_tbi) && !file.exists(idx_csi)) {
  message("Indexing VCF with bcftools...")
  status <- system(paste("bcftools index -t", shQuote(opt$vcf)))
  if (status != 0) stop("bcftools index failed")
}

gmmat <- paste0(opt$null_prefix, ".rda")
if (!file.exists(gmmat)) {
  # Newer SAIGE may write .rda or other suffixes; accept common patterns
  cands <- Sys.glob(paste0(opt$null_prefix, "*"))
  rda <- cands[grepl("\\.rda$", cands)]
  if (length(rda) < 1) stop(paste("Null model .rda not found for prefix", opt$null_prefix))
  gmmat <- rda[[1]]
}
var_ratio <- paste0(opt$null_prefix, ".varianceRatio.txt")
if (!file.exists(var_ratio)) {
  cands <- Sys.glob(paste0(opt$null_prefix, "*.varianceRatio.txt"))
  if (length(cands) < 1) stop("varianceRatio file not found")
  var_ratio <- cands[[1]]
}

# SAIGE >=0.38: --sampleFile is for BGEN only. For VCF, IDs come from the
# VCF header and are matched to the null model. If Step1 samples may be
# missing from the VCF, use --subSampleFile (optional).
subsample_arg <- ""
if (!is.na(opt$sample_file) && nzchar(opt$sample_file) && file.exists(opt$sample_file)) {
  subsample_arg <- paste0("--subSampleFile=", shQuote(opt$sample_file))
}

if (!file.exists(opt$sparse_grm)) stop(paste("sparse GRM not found:", opt$sparse_grm))
if (!file.exists(opt$sparse_grm_ids)) {
  stop(paste("sparse GRM sample IDs not found:", opt$sparse_grm_ids))
}

# Step2 has no --nThreads; pin BLAS/OpenMP to the WDL cpu allotment.
Sys.setenv(
  OMP_NUM_THREADS = as.character(opt$n_threads),
  OPENBLAS_NUM_THREADS = as.character(opt$n_threads),
  MKL_NUM_THREADS = as.character(opt$n_threads)
)
if (requireNamespace("RhpcBLASctl", quietly = TRUE)) {
  RhpcBLASctl::blas_set_num_threads(opt$n_threads)
  RhpcBLASctl::omp_set_num_threads(opt$n_threads)
}

cmd <- paste(
  "Rscript", shQuote(opt$step2_r),
  paste0("--vcfFile=", shQuote(opt$vcf)),
  paste0("--vcfFileIndex=", shQuote(if (file.exists(idx_tbi)) idx_tbi else idx_csi)),
  "--vcfField=GT",
  paste0("--chrom=", opt$chrom),
  paste0("--minMAC=", opt$min_mac),
  paste0("--GMMATmodelFile=", shQuote(gmmat)),
  paste0("--varianceRatioFile=", shQuote(var_ratio)),
  paste0("--SAIGEOutputFile=", shQuote(opt$out_raw)),
  paste0("--sparseGRMFile=", shQuote(opt$sparse_grm)),
  paste0("--sparseGRMSampleIDFile=", shQuote(opt$sparse_grm_ids)),
  subsample_arg,
  "--is_Firth_beta=TRUE",
  "--is_output_moreDetails=TRUE",
  "--LOCO=FALSE"
)
message("Running: ", cmd)
status <- system(cmd)
if (status != 0) stop(sprintf("SAIGE step2 failed with status %s", status))

raw <- fread(opt$out_raw, sep = "\t", header = TRUE, data.table = FALSE)
chr_c <- pick_col(raw, c("CHR", "Chromosome", "chrom", "chr"))
pos_c <- pick_col(raw, c("POS", "pos", "Position"))
id_c <- pick_col(raw, c("MarkerID", "SNPID", "ID", "SNP"))
a1_c <- pick_col(raw, c("Allele1", "A1", "REF"))
a2_c <- pick_col(raw, c("Allele2", "A2", "ALT"))
ac_c <- pick_col(raw, c("AC_Allele2", "AC", "MAC"))
af_c <- pick_col(raw, c("AF_Allele2", "AF", "MAF"))
beta_c <- pick_col(raw, c("BETA", "beta", "Beta"))
se_c <- pick_col(raw, c("SE", "se", "beta_SE"))
p_c <- pick_col(raw, c("p.value", "p.value.NA", "Pvalue", "P", "p"))
tstat_c <- pick_col(raw, c("Tstat", "Tstat_SPA", "SCORE"))
is_spa_c <- pick_col(raw, c("Is.SPA", "Is.SPA.converge", "SPA_converge"))
firth_c <- pick_col(raw, c("Is.FIRTH", "Is.FIRTH.converge", "FIRTH_converge"))

need <- c(chr_c, pos_c, beta_c, se_c, p_c)
if (any(is.na(need))) {
  stop(sprintf(
    "Could not map required SAIGE columns. Have: %s",
    paste(names(raw), collapse = ",")
  ))
}

n <- nrow(raw)
out <- data.frame(
  chrom = as.character(raw[[chr_c]]),
  pos = as.integer(raw[[pos_c]]),
  marker_id = if (!is.na(id_c)) as.character(raw[[id_c]]) else NA_character_,
  allele1 = if (!is.na(a1_c)) as.character(raw[[a1_c]]) else NA_character_,
  allele2 = if (!is.na(a2_c)) as.character(raw[[a2_c]]) else NA_character_,
  mac = if (!is.na(ac_c)) as.numeric(raw[[ac_c]]) else NA_real_,
  af = if (!is.na(af_c)) as.numeric(raw[[af_c]]) else NA_real_,
  beta = as.numeric(raw[[beta_c]]),
  se = as.numeric(raw[[se_c]]),
  pvalue = as.numeric(raw[[p_c]]),
  tstat = if (!is.na(tstat_c)) as.numeric(raw[[tstat_c]]) else NA_real_,
  is_spa = if (!is.na(is_spa_c)) as.character(raw[[is_spa_c]]) else NA_character_,
  is_firth = if (!is.na(firth_c)) as.character(raw[[firth_c]]) else NA_character_,
  stringsAsFactors = FALSE
)

# Prefer ALT/effect allele label for QC tables
out$allele <- ifelse(!is.na(out$allele2) & nzchar(out$allele2), out$allele2, out$allele1)

fwrite(out, opt$out_tsv, sep = "\t", na = "NA", quote = FALSE)
message(sprintf("Wrote normalized SAIGE results: %s (%d variants)", opt$out_tsv, nrow(out)))
