# Locityper (`locityper`)

Per-sample [Locityper](https://github.com/tprodanov/locityper) genotyping of
complex loci from a streamed CRAM subset. One workflow job = one sample.

Isaac’s Terra original: [EichlerLab/AoU_WDL `locityper_stream.wdl`](https://github.com/EichlerLab/AoU_WDL/blob/main/locityper/locityper_stream.wdl)
([`SOURCE.md`](SOURCE.md)). The WDL here is a **Verily Workbench** adaptation.
Do **not** submit it from Terra. In a VWB Jupyter app:

1. [`notebooks/rw/locityper_00_prep_reference.ipynb`](../notebooks/rw/locityper_00_prep_reference.ipynb)
   — NCBI GRCh38 no-alt FASTA, `samtools faidx`, Jellyfish 25-mers; locates the
   v9 CRAM manifest
2. [`notebooks/rw/locityper_01_run_stream.ipynb`](../notebooks/rw/locityper_01_run_stream.ipynb)
   — stage / register / submit via [`wb workflow`](https://support.workbench.verily.com/docs/guides/workflows/cromwell/)
   (same pattern as Matt’s `eTRs_getPhasedAlleleInfo` notebook)

The Workbench **Workflows GUI is unreliable**.

## What the WDL does

[`wdl/LocityperStream.wdl`](wdl/LocityperStream.wdl) — workflow name
`LocityperStream`.

| Step | What |
|---|---|
| `SplitBed` | Split the 4-column loci BED into shards of `N_split` lines |
| `LocityperPreprocessAndGenotype` (scatter) | Stream a windowed CRAM subset (`localization_optional` + `GCS_OAUTH_TOKEN`), `locityper preproc`, then `locityper genotype --subset-loci` per row via GNU parallel |
| `CombineTarFiles` | Merge per-shard `out_dir` tarballs |
| `Summarize` | Drop loci without `res.json.gz`; write a TSV of genotype / quality / read counts |

The genotype task also pulls Isaac’s GRCh38 background interval
(`chr17:72062001-76562000`) into the CRAM subset so Locityper can estimate
depth. That interval is **not** genotyped unless it is already in `bed`.
Override or clear `bg_region_bed` for a non-GRCh38 reference.

BED column 4 is the locus name (`cut -f4`). Missing names → empty
`--subset-loci` and failed shards.

Isaac’s Terra JSON used `--technology illumina` and a small CPU/memory
footprint (`locityper_n_cpu=2`). For PacBio HiFi set `technology` to whatever
`locityper preproc --help` lists (`hifi` / `pacbio`). Smoke-test **one
sample** and a tiny BED before a cohort CSV.

## Run on Verily Workbench

1. Open a JupyterLab app in the workspace (CLI `wb` is on `PATH`).
2. Clone this repo into the app if it is not already a git resource.
3. Run `locityper_00_prep_reference` (NCBI FASTA + Jellyfish). Size the VM
   to ≥16 GB RAM / ~20 GB free disk.
4. Open `locityper_01_run_stream`. Paste the `gs://` URIs and smoke-test
   CRAM from the prep notebook; add a loci BED + `vcf_db.tar.gz`.
5. Stage the WDL → register `locityper-stream` → submit. `SUBMIT` is off
   until you set it.

v9 Illumina CRAMs:
`workspace/vwb-aou-datasets-controlled-v9/v9/wgs/cram/manifest.csv`
(`gs://vwb-aou-datasets-controlled/pooled/wgs/cram/v8_base/wgs_{person_id}.cram`).

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
| `cram` / `crai` | Mapped WGS; CRAM is streamed, CRAI is localized |
| `ref_fa_uncompressed` / `ref_fai_uncompressed` | Uncompressed FASTA + fai (Isaac gunzipped assembly38 once and reused it) |
| `counts_jf` | Jellyfish k-mer counts for that reference (prep notebook: canonical 25-mers, `--lower-count 2`) |
| `bed` | `chrom start end locus_name` |
| `locityper_db_tar_gz` | `tar czf` of a Locityper DB directory named `vcf_db` (Isaac’s `locityper add -d vcf_db`) |
| `technology` | Default `illumina` |
| `locityper_docker` | Default `eichlerlab/locityper:1.4.5.0` (needs `gcloud`, `samtools` with libcurl, `parallel`, `locityper`) |
| `util_docker` | Default `python:3.11-slim` (split / tar / summarize) |

VPC-SC workspaces often cannot pull Docker Hub. Mirror both images to a
workspace-readable Artifact Registry and pass the AR tags.

## Outputs

| Output | What |
|---|---|
| `summary_csv` | `sample locus genotype quality total_reads unexpl_reads weight_dist warnings` |
| `results_tar_gz` | Combined `out_dir/loci/*/res.json.gz` |

Monitor with `wb workflow job list --workflow=locityper-stream` and
`wb workflow job describe --job-id=...`. Cancel with
`wb workflow job cancel --job-id=...`.
