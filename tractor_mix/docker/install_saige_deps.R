#!/usr/bin/env Rscript
# Lean SAIGE dependency installer for the pilot Docker image.
# Avoids SAIGE/extdata/install_packages.R (devtools/shiny stack) and GitHub
# API installs (rate-limited / moved repos on Cloud Build).

options(repos = c(CRAN = "https://cloud.r-project.org"))
options(Ncpus = max(1L, parallel::detectCores() - 1L))

need <- function(pkgs) {
  missing <- pkgs[!vapply(pkgs, requireNamespace, quietly = TRUE, FUN.VALUE = logical(1))]
  if (length(missing)) {
    install.packages(missing, dependencies = c("Depends", "Imports", "LinkingTo"))
  }
  for (p in pkgs) {
    if (!requireNamespace(p, quietly = TRUE)) {
      stop(sprintf("Failed to install package: %s", p), call. = FALSE)
    }
  }
}

# CRAN runtime / LinkingTo deps used by SAIGE 1.3.3
need(c(
  "remotes",
  "Rcpp",
  "RcppParallel",
  "RcppArmadillo",
  "RcppEigen",
  "data.table",
  "Matrix",
  "BH",
  "optparse",
  "RhpcBLASctl",
  "RSQLite",
  "dplyr",
  "SKAT",
  "R.utils",
  "MetaSKAT",   # now on CRAN (GitHub leeshawn/MetaSKAT is gone)
  "qlcMatrix",   # on CRAN; avoid GitHub API
  "slam",
  "sparsesvd",
  "docopt"
))

# SAIGE LinkingTo pins SPAtest == 3.1.2
if (!requireNamespace("SPAtest", quietly = TRUE) ||
    packageVersion("SPAtest") != "3.1.2") {
  remotes::install_version(
    "SPAtest",
    version = "3.1.2",
    repos = "https://cloud.r-project.org",
    upgrade = "never"
  )
}
if (packageVersion("SPAtest") != "3.1.2") {
  stop("SPAtest 3.1.2 is required", call. = FALSE)
}

message("SAIGE runtime dependencies OK")
print(sapply(
  c("Rcpp", "RcppArmadillo", "SPAtest", "SKAT", "MetaSKAT", "qlcMatrix", "RSQLite"),
  function(p) as.character(packageVersion(p))
))
