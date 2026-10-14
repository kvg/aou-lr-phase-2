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
# Exit status 1 if any phenotype is flagged.

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
              min_mu1mu = NA, coef_absmax = NA)
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
      if (length(why)) { rec$status <- "FLAGGED"; rec$reasons <- paste(why, collapse = "; ") }
    }
  }
  rows[[length(rows) + 1]] <- rec
  cat(sprintf("%-14s %-9s %s\n", rec$phenotype, rec$status, rec$reasons))
}

df <- do.call(rbind, lapply(rows, function(r) as.data.frame(r, stringsAsFactors = FALSE)))
write.table(df, out_path, sep = "\t", quote = FALSE, row.names = FALSE)

n_bad <- sum(df$status != "OK")
cat(sprintf("\n%d of %d null fits flagged. Wrote %s\n", n_bad, nrow(df), out_path))
if (n_bad > 0) {
  cat("Do NOT use a flagged phenotype's step-2 results: its score variance is wrong.\n")
  cat("Refit step 1 for those phenotypes (check covariate collinearity and\n")
  cat("case count) before rerunning the scan.\n")
}
quit(status = if (n_bad > 0) 1 else 0)
