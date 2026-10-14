#!/usr/bin/env Rscript
# Audit SAIGE/FELIX step-1 null fits for non-convergence before they reach step 2.
#
# Step 1 can return a diverged fit and still write a .rda, and nothing downstream
# notices: CheckNullArrays only compares array lengths, and the companion
# .varianceRatio.txt looks healthy either way. Observed on phenotype CA_142_2
# (phecode CA_142.2, 360 cases / 9,560 analysed): converged = FALSE,
# theta = (1, 2.04e13), intercept = -1.7e15, every fitted value underflowed to
# 2.22e-16 (= .Machine$double.eps). With mu*(1-mu) ~ 0 for every sample the score
# variance is ~0 at every marker, so the statistic diverges and step 2 dies with
#   Error in function boost::math::gamma_q_inv<double>(double, double): Overflow Error
#
# The crash is the LUCKY case. A less extreme non-converged fit does not abort --
# it returns p-values computed against a wrong variance, which is far worse
# because it is silent. Run this over every null before trusting a scan.
#
# Usage:
#   felix/scripts/audit_null_fits.R OUT.tsv NULL.rda [NULL.rda ...]
#   felix/scripts/audit_null_fits.R OUT.tsv $(ls nulls/*.null.rda)
#
# Flags a fit when any of:
#   converged is FALSE
#   theta has a component outside [1e-8, 1e6]
#   fitted values are constant, or outside (1e-10, 1-1e-10)
#   min mu*(1-mu) < 1e-8          (score variance collapses)
#   |coefficients| > 1e6
#   a low-cardinality covariate has a level with zero cases (separation)
# Status CHECK (not fatal) when the design is rank deficient but the fit is
# otherwise healthy: that is the expected shape for a correctly sex-restricted
# phenotype, whose global covariate list still carries a now-constant sex column.
# R's glm aliases such a column to NA and returns identical estimates either way,
# but SAIGE step 1 builds and inverts its own X'X, so it is reported rather than
# assumed safe. Exit status 1 if any phenotype is FLAGGED.
#
# The separation check is the one that names a CAUSE rather than a symptom. For
# CA_142.2 it reports sex: 4,295 participants of one sex hold 0 of the 360 cases,
# because the phecode is labelled "Both" in PhecodeX yet is female-only in
# practice, so selection never restricted it and sex stayed in the covariates.

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 2) {
  cat("usage: audit_null_fits.R OUT.tsv NULL.rda [NULL.rda ...]\n", file = stderr())
  quit(status = 2)
}
out_path <- args[1]
files <- args[-1]

THETA_LO <- 1e-8; THETA_HI <- 1e6
MU_EPS <- 1e-10; VAR_FLOOR <- 1e-8; COEF_HI <- 1e6

rows <- list()
for (f in files) {
  pheno <- sub("\\.null\\.rda$", "", basename(f))
  rec <- list(phenotype = pheno, file = f, status = "OK", reasons = "",
              converged = NA, n = NA, n_cases = NA,
              theta_max = NA, mu_min = NA, mu_max = NA, mu_distinct = NA,
              min_mu1mu = NA, coef_absmax = NA, separating_covariates = "",
              x_rank = NA, x_ncol = NA, constant_covariates = "")
  env <- new.env()
  loaded <- tryCatch({ load(f, envir = env); TRUE },
                     error = function(e) { rec$status <<- "UNREADABLE"
                                           rec$reasons <<- conditionMessage(e); FALSE })
  if (loaded) {
    if (!"modglmm" %in% ls(env)) {
      rec$status <- "UNREADABLE"; rec$reasons <- "no modglmm object"
    } else {
      m <- get("modglmm", envir = env)
      mu <- as.numeric(m$fitted.values)
      v  <- mu * (1 - mu)
      rec$converged   <- isTRUE(m$converged)
      rec$n           <- length(m$y)
      rec$n_cases     <- if (identical(m$traitType, "binary")) sum(m$y == 1) else NA
      rec$theta_max   <- max(abs(as.numeric(m$theta)))
      rec$mu_min      <- min(mu); rec$mu_max <- max(mu)
      rec$mu_distinct <- length(unique(mu))
      rec$min_mu1mu   <- min(v)
      rec$coef_absmax <- max(abs(as.numeric(m$coefficients)))

      why <- character(0)
      if (!isTRUE(m$converged)) why <- c(why, "converged=FALSE")
      th <- abs(as.numeric(m$theta))
      if (any(th > THETA_HI | th < THETA_LO)) why <- c(why, "theta out of range")
      if (rec$mu_distinct == 1)              why <- c(why, "fitted values constant")
      if (rec$mu_min <= MU_EPS || rec$mu_max >= 1 - MU_EPS) why <- c(why, "fitted values at 0/1")
      if (rec$min_mu1mu < VAR_FLOOR)         why <- c(why, "mu(1-mu) collapsed")
      if (rec$coef_absmax > COEF_HI)         why <- c(why, "coefficients diverged")

      # Design rank and constant columns.
      if (!is.null(m$X)) {
        Xm0 <- m$X
        rec$x_ncol <- ncol(Xm0)
        rec$x_rank <- tryCatch(qr(Xm0)$rank, error = function(e) NA_integer_)
        const <- character(0)
        for (j in seq_len(ncol(Xm0))) {
          if (length(unique(signif(Xm0[, j], 10))) == 1) {
            nm <- if (is.null(colnames(Xm0))) paste0("V", j) else colnames(Xm0)[j]
            if (!identical(nm, "minus1")) const <- c(const, nm)
          }
        }
        rec$constant_covariates <- paste(const, collapse = "; ")
      }

      # Separation: a covariate level holding zero cases. Only meaningful for
      # low-cardinality columns (sex, batch indicators), so continuous
      # covariates and the intercept are skipped.
      sep <- character(0)
      if (identical(m$traitType, "binary") && !is.null(m$X)) {
        Xm <- m$X; yv <- as.integer(m$y)
        for (j in seq_len(ncol(Xm))) {
          v <- Xm[, j]
          lv <- unique(signif(v, 10))
          if (length(lv) < 2 || length(lv) > 10) next
          key <- signif(v, 10)
          tab <- table(key, yv)
          if (!"1" %in% colnames(tab)) next
          zero_case <- rownames(tab)[tab[, "1"] == 0]
          if (length(zero_case)) {
            nm <- if (is.null(colnames(Xm))) paste0("V", j) else colnames(Xm)[j]
            n_zero <- sum(key %in% as.numeric(zero_case))
            sep <- c(sep, sprintf("%s(%d samples, 0 cases)", nm, n_zero))
          }
        }
      }
      if (length(sep)) {
        rec$separating_covariates <- paste(sep, collapse = "; ")
        why <- c(why, "covariate separates cases")
      }

      if (length(why)) {
        rec$status <- "FLAGGED"; rec$reasons <- paste(why, collapse = "; ")
      } else if (!is.na(rec$x_rank) && !is.na(rec$x_ncol) && rec$x_rank < rec$x_ncol) {
        rec$status <- "CHECK"
        rec$reasons <- sprintf("design rank %d of %d%s", rec$x_rank, rec$x_ncol,
                               if (nzchar(rec$constant_covariates))
                                 paste0(" (constant: ", rec$constant_covariates, ")") else "")
      }
    }
  }
  rows[[length(rows) + 1]] <- rec
  cat(sprintf("%-14s %-9s %s\n", rec$phenotype, rec$status, rec$reasons))
  if (nzchar(rec$separating_covariates))
    cat(sprintf("%-14s %-9s   separation: %s\n", "", "", rec$separating_covariates))
}

df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

n_bad <- sum(df$status == "FLAGGED" | df$status == "UNREADABLE")
n_chk <- sum(df$status == "CHECK")
cat(sprintf("\n%d of %d null fits flagged, %d to check. Wrote %s\n",
            n_bad, nrow(df), n_chk, out_path))
if (n_chk > 0) {
  cat("CHECK = rank-deficient design, fit otherwise healthy. Expected for a\n")
  cat("sex-restricted phenotype (constant sex column); confirm SAIGE step 1\n")
  cat("handled it rather than assuming so.\n")
}
if (n_bad > 0) {
  cat("Do NOT use a flagged phenotype's step-2 results: its score variance is wrong.\n")
  cat("Refit step 1 for those phenotypes before rerunning the scan. Where a\n")
  cat("covariate separates the cases, the fix is upstream in phenotype selection:\n")
  cat("restrict the phenotype to the affected sex/stratum and drop that covariate,\n")
  cat("rather than refitting with the same design.\n")
}
quit(status = if (n_bad > 0) 1 else 0)
