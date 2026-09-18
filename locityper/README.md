# Locityper (`locityper`)

Per-sample [Locityper](https://github.com/tprodanov/locityper) genotyping of
complex loci from AoU srWGS CRAMs. One workflow job = one sample. The WGS CRAM
stays in Nearline: a single `print_reads` pass builds a local minicram, then
`preproc` / `genotype` never touch GCS again.

Isaac’s Terra original: [EichlerLab/AoU_WDL `locityper_stream.wdl`](https://github.com/EichlerLab/AoU_WDL/blob/main/locityper/locityper_stream.wdl)
([`SOURCE.md`](SOURCE.md)). The WDL here is a **Verily Workbench** adaptation.
Do **not** submit it from Terra. In a VWB Jupyter app:

1. [`notebooks/rw/locityper_00_prep_reference.ipynb`](../notebooks/rw/locityper_00_prep_reference.ipynb)
   — GATK `Homo_sapiens_assembly38` FASTA + fai, Jellyfish 25-mers, one-locus
   smoke BED + toy `vcf_db.tar.gz`; locates the v9 CRAM manifest
2. [`notebooks/rw/locityper_01_run_stream.ipynb`](../notebooks/rw/locityper_01_run_stream.ipynb)
   — stage / register / submit via [`wb workflow`](https://support.workbench.verily.com/docs/guides/workflows/cromwell/)
   (same pattern as Matt’s `eTRs_getPhasedAlleleInfo` notebook)

The Workbench **Workflows GUI is unreliable**.

## What the WDL does

[`wdl/LocityperStream.wdl`](wdl/LocityperStream.wdl) — workflow name
`LocityperStream`.

| Step | What |
|---|---|
| `MakeMinicram` | Once per sample: [print_reads](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/print_reads.py) maps the loci BED + chr17 background through the CRAI and downloads each CRAM container **once** (GCS client, not `samtools view`) |
| `LocityperPreprocess` | Once per sample on that local minicram |
| `SplitBed` | Split the 4-column loci BED into shards of `N_split` lines |
| `LocityperGenotype` (scatter) | `locityper genotype --subset-loci` per row via GNU parallel; local minicram only |
| `CombineTarFiles` | Merge per-shard `out_dir` tarballs |
| `Summarize` | Drop loci without `res.json.gz`; write a TSV of genotype / quality / read counts |

`MakeMinicram` also pulls Isaac’s GRCh38 background interval
(`chr17:72062001-76562000`) so Locityper can estimate depth. That interval is
**not** genotyped unless it is already in `bed`. Override or clear
`bg_region_bed` for a non-GRCh38 reference. `minicram_n_preemptible` defaults
to **0** so a preempt does not re-read Nearline.

BED column 4 is the locus name (`cut -f4`). Missing names → empty
`--subset-loci` and failed shards.

Isaac’s Terra JSON used `--technology illumina` and a small CPU/memory
footprint (`locityper_n_cpu=2`). For PacBio HiFi set `technology` to whatever
`locityper preproc --help` lists (`hifi` / `pacbio`). Smoke-test **one
sample** and a tiny BED before a cohort CSV.

## Run on Verily Workbench

1. Open a JupyterLab app in the workspace (CLI `wb` is on `PATH`).
2. Clone this repo into the app if it is not already a git resource.
3. Build/push the print_reads image (`./build_docker.sh`) or mirror it into a
   workspace-readable Artifact Registry. VPC-SC cannot pull the default tag.
4. Run `locityper_00_prep_reference` (GATK `Homo_sapiens_assembly38` + Jellyfish
   + smoke catalog). Size the VM to ≥16 GB RAM / ~20 GB free disk.
5. Open `locityper_01_run_stream`. Paste the `gs://` URIs from the prep
   notebook (including `smoke.bed` + `vcf_db.smoke.tar.gz`). Set
   `PRINT_READS_DOCKER` if you mirrored the image. `RECREATE_WORKFLOW=True`
   once after this WDL change.
6. Stage the WDL → register `locityper-stream` → submit. `SUBMIT` is off
   until you set it.

v9 Illumina CRAMs:
`workspace/vwb-aou-datasets-controlled-v9/v9/wgs/cram/manifest.csv`
(`gs://vwb-aou-datasets-controlled/pooled/wgs/cram/v8_base/wgs_{person_id}.cram`).

Staged reference (prep notebook, 2026-09-15):

```
gs://aou-lr-phase2-resources/locityper/refs/Homo_sapiens_assembly38.fasta
gs://aou-lr-phase2-resources/locityper/refs/Homo_sapiens_assembly38.fasta.fai
gs://aou-lr-phase2-resources/locityper/refs/counts.Homo_sapiens_assembly38.k25.jf
gs://aou-lr-phase2-resources/locityper/smoke.bed
gs://aou-lr-phase2-resources/locityper/vcf_db.smoke.tar.gz
```

The smoke catalog is **not** a scientific haplotype panel. It is one CYP2D6
interval ([`configs/smoke.bed`](configs/smoke.bed), same span as
`bam_to_contig/regions/CYP2D6_through_7_GRCh38.bed`) with two haplotypes:
assembly38 plus a single interior SNP. Locityper refuses a one-haplotype DB,
so the SNP exists only to satisfy `locityper target`. Swap in Isaac’s real
`vcf_db.tar.gz` / loci BED before any analysis.

Single job:

```bash
wb workflow job run \
  --workflow=locityper-stream \
  --output-bucket-id="$OUTPUT_BUCKET_ID" \
  --output-path=workflowRuns/locityper \
  --inputs-uri=gs://BUCKET/locityper/inputs.json
```

Batch (one row per sample). Shared files must be repeated on every row;
Workbench does not fill WDL defaults from a sidecar JSON when you use a CSV:

```bash
wb workflow job run \
  --workflow=locityper-stream \
  --output-bucket-id="$OUTPUT_BUCKET_ID" \
  --output-path=workflowRuns/locityper \
  --batch-input-bucket-id="$OUTPUT_BUCKET_ID" \
  --batch-input-csv-path=locityper/batch.csv \
  --column-mapping-uri=gs://BUCKET/locityper/column_mapping.json
```

Examples: [`configs/locityper.inputs.json.example`](configs/locityper.inputs.json.example),
[`configs/batch.header.csv`](configs/batch.header.csv),
[`configs/column_mapping.json`](configs/column_mapping.json).

Register (after copying the WDL into a workspace bucket):

```bash
wb workflow create \
  --bucket-id="$OUTPUT_BUCKET_ID" \
  --path=wdl/locityper/LocityperStream.wdl \
  --workflow=locityper-stream \
  --workflow-type=WDL \
  --display-name="Locityper stream"
```

Git repos are **not** valid workflow sources; the WDL has to live in GCS.
See [Create batch jobs](https://support.workbench.verily.com/docs/guides/workflows/create_batch_jobs/).

## Inputs

| Input | Notes |
|---|---|
| `sample_id` | Output prefix |
| `cram` / `crai` | Mapped WGS. CRAM is `localization_optional` (stays in Nearline); CRAI is localized |
| `ref_fa_uncompressed` / `ref_fai_uncompressed` | Uncompressed FASTA + fai (AoU srWGS: GATK `Homo_sapiens_assembly38`) |
| `counts_jf` | Jellyfish k-mer counts for that reference (prep notebook: canonical 25-mers, `--lower-count 2`) |
| `bed` | `chrom start end locus_name`. Smoke: [`configs/smoke.bed`](configs/smoke.bed) (`CYP2D6_smoke`) |
| `locityper_db_tar_gz` | `tar czf` of a Locityper DB directory named `vcf_db`. Smoke: prep notebook `vcf_db.smoke.tar.gz` (toy 2-haplotype panel). Real runs need Isaac’s pangenome DB |
| `technology` | Default `illumina` |
| `gcloud_project` | Requester-pays billing project for the CRAM bucket; empty uses `$GOOGLE_CLOUD_PROJECT` |
| `print_reads_docker` | Default `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.0` ([`build_docker.sh`](build_docker.sh)) |
| `locityper_docker` | Default `eichlerlab/locityper:1.4.5.0` (needs `samtools`, `parallel`, `locityper`) |
| `util_docker` | Default `python:3.11-slim` (split / tar / summarize) |

VPC-SC workspaces often cannot pull Docker Hub or Broad AR. Mirror
`print_reads_docker` and `locityper_docker` (and `util_docker` if needed) to a
workspace-readable Artifact Registry.

## Outputs

| Output | What |
|---|---|
| `summary_csv` | `sample locus genotype quality total_reads unexpl_reads weight_dist warnings` |
| `results_tar_gz` | Combined `out_dir/loci/*/res.json.gz` |
| `minicram` / `minicrai` | Per-sample subset written by `print_reads` |
| `minicram_transfer_stats` | Bytes / containers pulled from Nearline |

Monitor with `wb workflow job list --workflow=locityper-stream` and
`wb workflow job describe --job-id=...`. Cancel with
`wb workflow job cancel --job-id=...`.
