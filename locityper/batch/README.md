# Locityper on native Google Batch

Same Nearline minicram → preproc/genotype path as
[`../wdl/LocityperStream.wdl`](../wdl/LocityperStream.wdl), submitted with
`gcloud batch` from a Verily Workbench Jupyter app. The WDL notebook
(`locityper_01_run_stream`) is unchanged.

dsub is **not** the path here (`submit_dsub.sh` remains as a leftover). Native
Batch uses a Python minicram worker in print-reads `0.1.2` (no bash/`gcloud`
in that image) and host `gcloud` I/O around the Locityper container.

The WGS CRAM stays a `gs://` env var (never localized). Minicram range-fetches
CRAM containers for the BED loci plus Isaac’s chr17 background interval.
Genotype localizes the **minicram** plus fasta / jellyfish / BED / `vcf_db`,
runs `locityper preproc` then per-locus `genotype` on one VM (no WDL scatter),
and uploads `gts.filtered.csv` + `locityper.tar.gz`.

VPC-SC cannot pull Docker Hub (`eichlerlab/locityper:1.4.5.0`). Set
`LOCITYPER_DOCKER` to a workspace-readable Artifact Registry mirror.

## Submit (VWB Jupyter)

`PET_SA_EMAIL` and `GOOGLE_CLOUD_PROJECT` are set by the app. Prefer
[`../../notebooks/rw/locityper_02_run_batch.ipynb`](../../notebooks/rw/locityper_02_run_batch.ipynb)
(`SUBMIT_*` off until you flip one). Default BED / `vcf_db` are the toy CYP2D6
smoke catalog from `locityper_00_prep_reference`. Swap those URIs before any
analysis cohort.

Cohort runs use a GCS ledger (`scripts/vwb_batch/`):
[`locityper_03_budget.ipynb`](../../notebooks/rw/locityper_03_budget.ipynb)
then
[`locityper_04_dispatch.ipynb`](../../notebooks/rw/locityper_04_dispatch.ipynb).
Read-only status:
[`batch_monitor.ipynb`](../../notebooks/rw/batch_monitor.ipynb) (`JOB_NAME_PREFIX = "lt-"`).
Tasks read a GCS TSV (`TASKS_TSV` + `BATCH_TASK_INDEX`) so job JSON stays
under 1 MiB. Pilot workers write `*.resources.tsv`; `_03` tunes compute into
`run.json` for `_04`. Or:

```bash
./submit_batch.sh --csv ../configs/batch.header.csv --stage minicram
./submit_batch.sh --csv ../configs/batch.header.csv --stage genotype
```

## Outputs

Under `--out-prefix` (default `$OUTPUT_BUCKET_GS/batchRuns/locityper`):

```
{sample_id}/{sample_id}.minicram.cram
{sample_id}/{sample_id}.minicram.cram.crai
{sample_id}/{sample_id}.data_transfer_stats.tsv
{sample_id}/{sample_id}.gts.filtered.csv
{sample_id}/{sample_id}.locityper.tar.gz
```

Run ledger: `batchRuns/locityper/runs/{RUN_ID}/`.

The Cromwell notebooks (`locityper_01_run_stream`) are unchanged.
