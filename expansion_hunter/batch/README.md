# ExpansionHunter on native Google Batch

Same minicram → genotype split as
[`../wdl/ExpansionHunterMinicram.wdl`](../wdl/ExpansionHunterMinicram.wdl),
submitted with `gcloud batch` from a Verily Workbench Jupyter app. The WDL
notebook (`expansion_hunter_01_run`) is unchanged.

dsub is **not** the path here. Native Batch works: host GCS smoke, in-container
`google.cloud.storage` (`eh-ctr-smoke-260929-034427`), two-locus minicram
(`eh-mini-260929-035922`), and the 711-locus degenerate panel on sample
`1000000` (`eh-mini-260929-051752` ~8 min / 67.5 MiB, `eh-gt-260929-053643`
~1 min).

The WGS CRAM stays a `gs://` env var (never localized). Minicram downloads
catalog, hg38, and the CRAI with the GCS client, then range-fetches CRAM
containers. Genotype localizes the **minicram** (small) plus fasta/catalog
with host `gcloud storage`, runs ExpansionHunter in the EH image (no GCS
client in `0.1.0`), and uploads JSON/VCF. Do not call `gcloud` inside
print-reads `0.1.2` (`gcloud`’s shebang is bash; slim returns exit 127).

## Submit (VWB Jupyter)

`PET_SA_EMAIL` and `GOOGLE_CLOUD_PROJECT` are set by the app. Prefer
[`../../notebooks/rw/expansion_hunter_02_run_batch.ipynb`](../../notebooks/rw/expansion_hunter_02_run_batch.ipynb)
(`SUBMIT_*` off until you flip one). `CATALOG_KIND = "degenerate"` is the
711-locus production panel (bigger smoke on sample `1000000`). `"tiny"` is
the two-locus path.

Cohort runs use a GCS ledger (`scripts/vwb_batch/`):
[`expansion_hunter_03_budget.ipynb`](../../notebooks/rw/expansion_hunter_03_budget.ipynb)
then
[`expansion_hunter_04_dispatch.ipynb`](../../notebooks/rw/expansion_hunter_04_dispatch.ipynb).
Read-only status of smoke + ledger jobs:
[`batch_monitor.ipynb`](../../notebooks/rw/batch_monitor.ipynb).
Tasks read a GCS TSV (`TASKS_TSV` + `BATCH_TASK_INDEX`) so job JSON stays
under 1 MiB at 100k samples. Pilot workers write peak RSS/CPU to
`*.resources.tsv`; `_03` tunes `cpuMilli` / `memoryMib` / `bootDiskMib` into
`run.json` for `_04`. Slide figures from the 711-locus smoke live in
[`figures/`](figures/). Or:

```bash
./submit_batch.sh --csv ../configs/batch.header.csv --stage minicram
# after minicram SUCCEEDED:
./submit_batch.sh --csv ../configs/batch.header.csv --stage genotype
```

Watch:

```bash
gcloud batch jobs describe JOB_ID --location=us-central1 \
  --format='yaml(status.state,status.statusEvents)'
```

On SUCCESS:

```
gs://…/batchRuns/expansion_hunter/{sample_id}/{sample_id}.minicram.cram
gs://…/batchRuns/expansion_hunter/{sample_id}/{sample_id}.minicram.cram.crai
gs://…/batchRuns/expansion_hunter/{sample_id}/{sample_id}.data_transfer_stats.tsv
gs://…/batchRuns/expansion_hunter/{sample_id}/{sample_id}.minicram.worker.log
gs://…/batchRuns/expansion_hunter/{sample_id}/{sample_id}.EH.json
gs://…/batchRuns/expansion_hunter/{sample_id}/{sample_id}.EH.vcf
gs://…/batchRuns/expansion_hunter/{sample_id}/{sample_id}.EH.worker.log
```

PET cannot read Cloud Logging. Use the worker log objects.

Cancel: `gcloud batch jobs delete JOB_ID --location=us-central1`.
The old [`submit_dsub.sh`](submit_dsub.sh) remains for reference only.
