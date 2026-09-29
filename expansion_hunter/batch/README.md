# ExpansionHunter on Google Cloud Batch (dsub)

Same minicram → genotype split as
[`../wdl/ExpansionHunterMinicram.wdl`](../wdl/ExpansionHunterMinicram.wdl),
submitted with [`dsub --provider google-batch`](https://github.com/DataBiosphere/dsub)
instead of Cromwell. Use this from a Verily Workbench Jupyter app when
Cromwell is flaky. VWB’s dsub flags (PET SA + private VPC) are documented
at [Get started with dsub on Verily Workbench](https://support.workbench.verily.com/docs/guides/workflows/dsub/).

The WGS CRAM is an `--env` `gs://` URI (never `--input`), so Batch does not
localize the Nearline object. The CRAI is `--env` `gs://` as well (same
requester-pays bucket). Stage 1 (`make_minicram.sh`) runs in the
print_reads image; stage 2 (`genotype.sh`) runs in the EH image. Two jobs
because the images differ; `--after` chains them.

## Submit (VWB Jupyter)

`dsub` must be on `PATH` (`pip install dsub` if it is not). `PET_SA_EMAIL`
and `GOOGLE_CLOUD_PROJECT` are set by the app. VPC-SC still needs
workspace-readable image mirrors.

```bash
# one-row smoke CSV: configs/batch.header.csv
./submit_dsub.sh --csv ../configs/batch.header.csv --stage minicram

# after minicram succeeds:
./submit_dsub.sh --csv ../configs/batch.header.csv --stage genotype \
  --after MINICRAM_JOB_ID

# smoke only: --stage all waits on minicram via --after, then submits genotype
./submit_dsub.sh --csv ../configs/batch.header.csv --stage all
```

Do **not** `--stage all` on a cohort CSV from a notebook cell unless you
want the kernel blocked until every minicram task finishes.

Notebook:
[`../../notebooks/rw/expansion_hunter_02_run_batch.ipynb`](../../notebooks/rw/expansion_hunter_02_run_batch.ipynb)
(`SUBMIT` off until you flip it). Prep refs with
`expansion_hunter_00_prep_reference` first.

## Outputs

Under `--out-prefix` (default `$OUTPUT_BUCKET_GS/batchRuns/expansion_hunter`):

```
{sample_id}/{sample_id}.minicram.cram
{sample_id}/{sample_id}.minicram.cram.crai
{sample_id}/{sample_id}.data_transfer_stats.tsv
{sample_id}/{sample_id}.EH.json
{sample_id}/{sample_id}.EH.vcf
logs/
```

## Status

```bash
dstat --provider google-batch --project "$GOOGLE_CLOUD_PROJECT" \
  --jobs JOB_ID --status '*'
```

Cancel: `ddel --provider google-batch --project "$GOOGLE_CLOUD_PROJECT" --jobs JOB_ID`.
The Cromwell notebooks (`expansion_hunter_01_run`) are unchanged.
