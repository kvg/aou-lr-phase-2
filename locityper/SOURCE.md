Vendored from [EichlerLab/AoU_WDL `locityper/locityper_stream.wdl`](https://github.com/EichlerLab/AoU_WDL/blob/main/locityper/locityper_stream.wdl)
at `7c7bf61151f501cded419d5e1d74bc8726298947` (Isaac Wong, 2026-07-24, “reduce log outputs from locityper”).

The WDL in this repo (`wdl/LocityperStream.wdl`) is a **Verily Workbench**
adaptation of that per-sample stream/genotype workflow. Differences from
upstream:

- Workflow renamed `ValidateVariants` → `LocityperStream` (input JSON keys
  change).
- Dropped unused helper tasks (`GenerateDBFromVCF`, `GunzipReference`,
  `SplitBedNames`, `SubsetBed`, `FilterNames`, `FilterBed`).
- Dropped the `docker://` image prefix (GCP Cromwell).
- `locityper_docker` / `print_reads_docker` / `util_docker` are workflow
  inputs so VPC-SC workspaces can point at Artifact Registry mirrors.
- `technology` and `bg_region_bed` are inputs (Isaac hardcoded
  `--technology illumina` and the GRCh38 chr17 background interval).
- `Summarize` inlines JSON→TSV instead of `/locityper/extra/into_csv.py` on
  the Broad `lr-hidive` image (that GCR tag is often not pullable from VWB).
- CRAM extract is **not** `samtools view` over `GCS_OAUTH_TOKEN`. Isaac’s
  combined `LocityperPreprocessAndGenotype` scatter streamed the Nearline CRAM
  (and the 4.5 Mb chr17 background) once per shard. Here `MakeMinicram` runs
  [str-analysis `print_reads`](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/print_reads.py)
  once per sample (CRAI → unique containers), `LocityperPreprocess` runs once
  on that minicram, and `LocityperGenotype` scatters without GCS.

Isaac’s Terra input JSON (data-table placeholders, `fc-secure-` URIs) lived at
[`locityper/inputs.locityper.json`](https://github.com/EichlerLab/AoU_WDL/blob/main/locityper/inputs.locityper.json).
The Dockerfile next to it builds `eichlerlab/locityper:1.4.5.0` (gcloud +
samtools `--enable-libcurl` + GNU parallel).
