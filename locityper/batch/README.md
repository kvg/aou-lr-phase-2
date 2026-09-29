# Locityper on Google Cloud Batch (dsub)

Same Nearline minicram → preproc/genotype path as
[`../wdl/LocityperStream.wdl`](../wdl/LocityperStream.wdl), submitted with
[`dsub --provider google-batch`](https://github.com/DataBiosphere/dsub)
instead of Cromwell. Use this from a Verily Workbench Jupyter app when
Cromwell is flaky. VWB’s dsub flags (PET SA + private VPC) are documented
at [Get started with dsub on Verily Workbench](https://support.workbench.verily.com/docs/guides/workflows/dsub/).

The WGS CRAM is an `--env` `gs://` URI (never `--input`). The CRAI is `--env`
`gs://` as well (same requester-pays bucket). Stage 1
(`make_minicram.sh`) is `print_reads` plus Isaac’s chr17 background interval.
Stage 2 (`genotype.sh`) runs `locityper preproc`, then GNU parallel over
every BED locus on **one** VM (no WDL scatter), then the same summary TSV.

## Submit (VWB Jupyter)

On a VWB Jupyter terminal, run `dsub_activate` first ([VWB dsub
guide](https://support.workbench.verily.com/docs/guides/workflows/dsub/)).
dsub google-batch copies `--input`/`--logging` with `DSUB_CLOUD_SDK_IMAGE`
(see [`vwb_dsub.sh`](../../scripts/vwb_dsub.sh)); VPC-SC cannot pull the
default `gcr.io/.../cloud-sdk:slim`.

```bash
./submit_dsub.sh --csv ../configs/batch.header.csv --stage minicram
./submit_dsub.sh --csv ../configs/batch.header.csv --stage genotype \
  --after MINICRAM_JOB_ID
```

`--stage all` waits on minicram via `--after` before submitting genotype —
fine for the one-row smoke CSV, not for a cohort from a notebook cell.

Notebook:
[`../../notebooks/rw/locityper_02_run_batch.ipynb`](../../notebooks/rw/locityper_02_run_batch.ipynb)
(`SUBMIT` off). Prep refs with `locityper_00_prep_reference` first.

`LOCITYPER_DOCKER` defaults to `eichlerlab/locityper:1.4.5.0`. VPC-SC cannot
pull Docker Hub; mirror it (and the print_reads image) into a
workspace-readable Artifact Registry and set `PRINT_READS_DOCKER` /
`LOCITYPER_DOCKER`.

## Outputs

Under `--out-prefix` (default `$OUTPUT_BUCKET_GS/batchRuns/locityper`):

```
{sample_id}/{sample_id}.minicram.cram
{sample_id}/{sample_id}.minicram.cram.crai
{sample_id}/{sample_id}.data_transfer_stats.tsv
{sample_id}/{sample_id}.gts.filtered.csv
{sample_id}/{sample_id}.locityper.tar.gz
logs/
```

## Status

```bash
dstat --provider google-batch --project "$GOOGLE_CLOUD_PROJECT" \
  --location us-central1 --users jupyter --jobs JOB_ID --status '*'
```

dsub `--logging` objects are `{job-id}.{task-id}.log` (plus `-stdout` /
`-stderr`) under `logs/`, not `{sample}.minicram.worker.log`. If those are
missing:

```bash
gcloud batch jobs describe JOB_ID-1-0 \
  --project "$GOOGLE_CLOUD_PROJECT" --location us-central1
```

The Cromwell notebooks (`locityper_01_run_stream`) are unchanged.
