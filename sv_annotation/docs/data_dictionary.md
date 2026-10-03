# Data dictionary

## Unified site table columns

| Column | Meaning |
|--------|---------|
| chrom, pos, end | GRCh38 coordinates; pos is 1-based VCF POS; end from INFO/END or length |
| id | Variant ID (synthesized if missing) |
| ref, alt | VCF alleles; long sequence ALTs are stored as `<SVTYPE>` (CADD-SV uses coordinates, not allele sequence) |
| svtype | Canonical class (DEL/DUP/INS/INV/BND/…) |
| svlen | Signed length when available; empty for unresolved BNDs |
| filter | VCF FILTER |
| source_vcf | `main`, `bnd`, or `large` |
| phase | `phase1` or `phase2` |
| ac, an, af | Cohort allele count/number/frequency; no-calls excluded |
| n_carriers | Samples with ≥1 alt allele called |
| freq_class | singleton / polymorphic / major / shared |
| region_class | Exclusive context class: `US`, `RM`, `SD`, or `SR` (SR overrides SD overrides RM). Breakpoint rule; body coverage >50% for DEL/DUP/CNV longer than 5 kb |
| hit_* | Evidence before the priority collapse (`hit_rmsk`, `hit_simpleRepeat`, `hit_genomicSuperDups` can co-occur). `hit_cmrg` is reference-span overlap and is not part of `region_class` |
| cadd_sv_phred | CADD-SV PHRED (empty if unscored) |
| cadd_sv_bin | low / mid / high / unscored (`large` partition and DUP/INV > 1 Mb are left unscored) |
| size_bin_ge20, size_bin_ge50 | Absolute length thresholds |
| suppressed_main_duplicate | true only on rows dropped from main due to large overlap |

## Manuscript counts TSV

| Column | Meaning |
|--------|---------|
| metric | deletions / deletions_{us,rm,sd,sr,cmrg} / deletions_context_pct_{ge50,ge20} / duplications / insertions / insertions_{us,rm,sd,sr,cmrg} / insertions_context_pct_{ge50,ge20} / inversions / breakends / large_events_gt10kb. Parent class row is ALL. `*_context_pct_ge50` and `*_context_pct_ge20` are US/RM/SD/SR/CMRG as percents of sites at that size cutoff. US+RM+SD+SR = 100 within the row; CMRG is independent. |
| phase*_display | `ge50; ge20` for resolved classes, or `—` if unavailable |
| phase*_ge50 / ge20 | Numeric counts (empty when unavailable) |

## Discovery TSV

| Column | Meaning |
|--------|---------|
| sample_index | 0-based position in ancestry-ordered sample list |
| sample_id | Sample identifier |
| ancestry | Ancestry label used for ordering |
| stratum | Stratification key (`all`, region class, and/or CADD bin) |
| freq_class | Frequency stack layer |
| cumulative_count | Sites of that class discovered by this sample index |
