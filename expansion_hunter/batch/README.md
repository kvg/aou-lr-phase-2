# ExpansionHunter on native Google Batch

Same minicram → genotype split as
[`../wdl/ExpansionHunterMinicram.wdl`](../wdl/ExpansionHunterMinicram.wdl),
submitted with `gcloud batch` from a Verily Workbench Jupyter app. The WDL
notebook (`expansion_hunter_01_run`) is unchanged.

dsub is **not** the path here. It is Batch plus sidecars that cannot pull
`gcr.io` cloud-sdk inside AoU VPC-SC, so logs and `--input` copies never
appear. Native Batch with one print-reads container works: host GCS smoke
and in-container `google.cloud.storage` both SUCCEEDED
(`eh-ctr-smoke-260929-034427`).

The WGS CRAM stays a `gs://` env var (never localized). The worker
downloads catalog, hg38, and the CRAI with the GCS client, then
`make_minicram_for_expansion_hunter` range-fetches CRAM containers. Do not
use `gcloud` inside this image until a later tag installs bash (`gcloud`’s
shebang is bash; slim 0.1.2 returns exit 127).

## Submit (VWB Jupyter)

`PET_SA_EMAIL` and `GOOGLE_CLOUD_PROJECT` are set by the app.

```bash
./submit_batch.sh --csv ../configs/batch.header.csv --stage minicram
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
```

PET cannot read Cloud Logging. Use the worker log object, not
`gcloud logging`.

Genotype (EH image) is not wired yet — submit it after the minicram objects
exist. The old [`submit_dsub.sh`](submit_dsub.sh) remains for reference only.

Cancel: `gcloud batch jobs delete JOB_ID --location=us-central1`.
