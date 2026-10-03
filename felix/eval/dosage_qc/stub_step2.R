# Stub for step2_SPAtests.R: records args + FELIX_DOSAGE_QC and writes a minimal result.
args <- commandArgs(trailingOnly = TRUE)
writeLines(c(paste("ENV", Sys.getenv("FELIX_DOSAGE_QC", "<unset>")), args), "val/stub_args.txt")
out <- sub("^--SAIGEOutputFile=", "", grep("^--SAIGEOutputFile=", args, value = TRUE))
out <- gsub("^'|'$", "", out)
writeLines(c("CHR\tPOS\tP_cct_admixed_c\tP_het_admixed_c\tP_hom_admixed_c", "chr1\t1\t0.5\t0.5\t0.5"), out)
