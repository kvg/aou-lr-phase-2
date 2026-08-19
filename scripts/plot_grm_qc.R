#!/usr/bin/env Rscript
# QC plots / tables for a Tractor-Mix sparse GRM (+ optional pedigree covariates).
#
# Typical Workbench usage after a successful MakeGRM:
#   Rscript plot_grm_qc.R \
#     --sparse-rds grm_sparse.rds \
#     --covariates covariates.source_rebuilt.csv.gz \
#     --out-dir grm_qc

suppressPackageStartupMessages({
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
parse_args <- function(args) {
  out <- list(
    sparse_rds = NA_character_,
    bands = NA_character_,
    histogram = NA_character_,
    close_pairs = NA_character_,
    covariates = NA_character_,
    threshold = 0.05,
    out_dir = "grm_qc"
  )
  i <- 1
  while (i <= length(args)) {
    key <- args[[i]]
    known <- c(
      "--sparse-rds", "--bands", "--histogram", "--close-pairs",
      "--covariates", "--threshold", "--out-dir"
    )
    if (key %in% known) {
      if (i == length(args)) stop(paste("Missing value for", key))
      val <- args[[i + 1]]
      if (key == "--sparse-rds") out$sparse_rds <- val
      if (key == "--bands") out$bands <- val
      if (key == "--histogram") out$histogram <- val
      if (key == "--close-pairs") out$close_pairs <- val
      if (key == "--covariates") out$covariates <- val
      if (key == "--threshold") out$threshold <- as.numeric(val)
      if (key == "--out-dir") out$out_dir <- val
      i <- i + 2
    } else {
      stop(paste("Unknown arg:", key))
    }
  }
  out
}

opt <- parse_args(args)
dir.create(opt$out_dir, showWarnings = FALSE, recursive = TRUE)

# --- Prefer precomputed TSVs from sparsify_grm.R; else derive from sparse RDS ---
hist_df <- NULL
bands_df <- NULL
close_df <- NULL
K <- NULL
grm_ids <- NULL

if (!is.na(opt$histogram) && file.exists(opt$histogram)) {
  hist_df <- read.delim(opt$histogram, stringsAsFactors = FALSE)
}
if (!is.na(opt$bands) && file.exists(opt$bands)) {
  bands_df <- read.delim(opt$bands, stringsAsFactors = FALSE)
}
if (!is.na(opt$close_pairs) && file.exists(opt$close_pairs)) {
  close_df <- read.delim(opt$close_pairs, stringsAsFactors = FALSE)
}

if (!is.na(opt$sparse_rds) && file.exists(opt$sparse_rds)) {
  message("Reading sparse RDS: ", opt$sparse_rds)
  K <- readRDS(opt$sparse_rds)
  grm_ids <- rownames(K)
  if (is.null(grm_ids)) grm_ids <- as.character(seq_len(nrow(K)))
}

if ((is.null(hist_df) || is.null(bands_df) || is.null(close_df)) && !is.null(K)) {
  message("Deriving QC tables from sparse RDS")
  ids <- grm_ids
  tri <- summary(Matrix::triu(K, k = 1))
  tri <- tri[tri$x >= opt$threshold, , drop = FALSE]
  kin <- tri$x
  if (is.null(hist_df) && length(kin) > 0) {
    kin_max <- max(kin, na.rm = TRUE)
    br <- seq(opt$threshold, kin_max, by = 0.025)
    if (length(br) < 2 || tail(br, 1) <= kin_max) {
      br <- c(br, kin_max + 0.025)
    }
    br <- sort(unique(br))
    h <- hist(kin, breaks = br, plot = FALSE, right = FALSE, include.lowest = TRUE)
    hist_df <- data.frame(
      kinship_bin_lo = h$breaks[-length(h$breaks)],
      kinship_bin_hi = h$breaks[-1],
      n_pairs = as.integer(h$counts)
    )
  }
  if (is.null(bands_df) && length(kin) > 0) {
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
    bands_df <- as.data.frame(table(band = b), stringsAsFactors = FALSE)
    names(bands_df)[2] <- "n_pairs"
    bands_df$expected_relationship <- expected[match(bands_df$band, labels)]
    bands_df$frac_of_pairs_ge_0.05 <- bands_df$n_pairs / length(kin)
  }
  if (is.null(close_df)) {
    keep <- tri$x >= 0.2
    close_df <- data.frame(
      id1 = ids[tri$i[keep]],
      id2 = ids[tri$j[keep]],
      kinship = tri$x[keep],
      stringsAsFactors = FALSE
    )
    close_df <- close_df[order(-close_df$kinship), , drop = FALSE]
  }
}

if (is.null(hist_df) || is.null(bands_df)) {
  stop("Need --histogram/--bands from MakeGRM, or a readable --sparse-rds")
}

write.table(
  bands_df,
  file.path(opt$out_dir, "grm_relationship_bands.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
write.table(
  hist_df,
  file.path(opt$out_dir, "grm_kinship_histogram.tsv"),
  sep = "\t", quote = FALSE, row.names = FALSE
)
if (!is.null(close_df)) {
  write.table(
    close_df,
    file.path(opt$out_dir, "grm_close_pairs.tsv"),
    sep = "\t", quote = FALSE, row.names = FALSE
  )
}

# --- Pedigree ↔ GRM (precision + recall) ---
if (!is.na(opt$covariates) && file.exists(opt$covariates) && !is.null(K)) {
  message("Pedigree ↔ GRM checks using: ", opt$covariates)
  cov <- read.csv(opt$covariates, stringsAsFactors = FALSE, check.names = FALSE)
  id_col <- if ("research_id" %in% names(cov)) "research_id" else names(cov)[[1]]
  cov[[id_col]] <- as.character(cov[[id_col]])
  if (!("pedigree_family_id" %in% names(cov))) {
    stop("covariates lack pedigree_family_id")
  }
  cov$pedigree_family_id <- as.character(cov$pedigree_family_id)
  cov$pedigree_family_id[is.na(cov$pedigree_family_id) | cov$pedigree_family_id == ""] <-
    NA_character_
  fam_map <- setNames(cov$pedigree_family_id, cov[[id_col]])

  # A) Genetic close pairs -> same pedigree family (precision-like)
  if (!is.null(close_df) && nrow(close_df) > 0) {
    close_df$fam1 <- unname(fam_map[as.character(close_df$id1)])
    close_df$fam2 <- unname(fam_map[as.character(close_df$id2)])
    close_df$same_pedigree_family <-
      !is.na(close_df$fam1) & !is.na(close_df$fam2) & (close_df$fam1 == close_df$fam2)
    write.table(
      close_df,
      file.path(opt$out_dir, "grm_close_pairs_with_pedigree.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    ped_precision <- data.frame(
      direction = "genetic_close_to_pedigree",
      n_pairs = nrow(close_df),
      n_same_pedigree_family = sum(close_df$same_pedigree_family, na.rm = TRUE),
      frac_same_pedigree_family = mean(close_df$same_pedigree_family, na.rm = TRUE),
      stringsAsFactors = FALSE
    )
    message("Among GRM close pairs (kinship >= 0.2), share pedigree_family_id:")
    print(ped_precision)
    write.table(
      ped_precision,
      file.path(opt$out_dir, "grm_pedigree_concordance.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
  }

  # B) Recorded same-family pairs in GRM cohort -> kinship (recall)
  cov_in_grm <- data.frame(
    id = cov[[id_col]],
    fam = cov$pedigree_family_id,
    stringsAsFactors = FALSE
  )
  cov_in_grm <- cov_in_grm[cov_in_grm$id %in% grm_ids & !is.na(cov_in_grm$fam), , drop = FALSE]
  fam_sizes <- as.data.frame(table(fam = cov_in_grm$fam), stringsAsFactors = FALSE)
  names(fam_sizes)[2] <- "n_in_grm"
  multi <- fam_sizes$fam[fam_sizes$n_in_grm >= 2]
  message(sprintf(
    "People with pedigree_family_id in GRM: %d; families with >=2 in GRM: %d",
    nrow(cov_in_grm), length(multi)
  ))

  if (length(multi) > 0) {
    ped_pairs <- lapply(multi, function(fam) {
      ids <- cov_in_grm$id[cov_in_grm$fam == fam]
      comb <- utils::combn(ids, 2)
      data.frame(
        id1 = comb[1, ],
        id2 = comb[2, ],
        pedigree_family_id = fam,
        stringsAsFactors = FALSE
      )
    })
    ped_pairs_df <- do.call(rbind, ped_pairs)
    # Vectorized sparse lookup; 0 => below sparsify threshold
    ped_pairs_df$kinship <- as.numeric(K[cbind(ped_pairs_df$id1, ped_pairs_df$id2)])
    ped_pairs_df$below_sparse_threshold <- ped_pairs_df$kinship < opt$threshold
    ped_pairs_df <- ped_pairs_df[order(-ped_pairs_df$kinship), , drop = FALSE]
    write.table(
      ped_pairs_df,
      file.path(opt$out_dir, "grm_pedigree_pairs_kinship.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )

    thresholds <- c(0.05, 0.088, 0.177, 0.2, 0.354, 0.707)
    ped_recall <- do.call(rbind, lapply(thresholds, function(th) {
      data.frame(
        direction = "pedigree_to_grm",
        kinship_threshold = th,
        n_pedigree_pairs_in_grm = nrow(ped_pairs_df),
        n_families_ge2_in_grm = length(multi),
        n_recovered = sum(ped_pairs_df$kinship >= th),
        frac_recovered = mean(ped_pairs_df$kinship >= th),
        median_kinship = stats::median(ped_pairs_df$kinship),
        mean_kinship = mean(ped_pairs_df$kinship),
        stringsAsFactors = FALSE
      )
    }))
    write.table(
      ped_recall,
      file.path(opt$out_dir, "grm_pedigree_recall.tsv"),
      sep = "\t", quote = FALSE, row.names = FALSE
    )
    message("Recorded pedigree pairs recovered by GRM (this answers recapitulation):")
    print(ped_recall)
  } else {
    message("No pedigree families with >=2 members overlap the GRM cohort.")
  }
} else if (!is.na(opt$covariates) && file.exists(opt$covariates) && is.null(K)) {
  message("NOTE: need --sparse-rds for pedigree-pair recall")
}

# --- Plots ---
png(file.path(opt$out_dir, "grm_kinship_histogram.png"), width = 900, height = 600)
mid <- (hist_df$kinship_bin_lo + hist_df$kinship_bin_hi) / 2
bp <- barplot(
  hist_df$n_pairs,
  names.arg = sprintf("%.2f", mid),
  las = 2,
  cex.names = 0.7,
  main = "Pairwise kinship (pairs >= sparsify threshold)",
  xlab = "Kinship bin midpoint (PLINK --make-rel / GCTA-like scale)",
  ylab = "Number of pairs",
  col = "#4C78A8"
)
abline(v = approx(mid, bp, xout = c(0.088, 0.177, 0.354, 0.707))$y, lty = 2, col = "gray40")
legend(
  "topright",
  legend = c("~3rd", "~2nd", "~1st", "dup/MZ"),
  lty = 2,
  col = "gray40",
  bty = "n",
  title = "Band edges"
)
dev.off()

png(file.path(opt$out_dir, "grm_relationship_bands.png"), width = 900, height = 600)
par(mar = c(10, 5, 4, 2))
barplot(
  bands_df$n_pairs,
  names.arg = paste0(bands_df$band, "\n", bands_df$expected_relationship),
  las = 2,
  cex.names = 0.65,
  main = "Relatedness bands among pairs with kinship >= 0.05",
  ylab = "Number of pairs",
  col = c("#9ECAE1", "#6BAED6", "#4292C6", "#2171B5", "#08306B")
)
dev.off()

message("Wrote figures under ", normalizePath(opt$out_dir))
message("Bands:")
print(bands_df)
quit(save = "no", status = 0)
