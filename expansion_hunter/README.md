# ExpansionHunter (`expansion_hunter`)

Per-sample [ExpansionHunter](https://github.com/bw2/ExpansionHunter) genotyping
of tandem-repeat loci from AoU srWGS CRAMs. One workflow job = one sample. The
WGS CRAM stays in Nearline: a single
[`make_minicram_for_expansion_hunter`](https://github.com/broadinstitute/str-analysis/blob/main/str_analysis/make_minicram_for_expansion_hunter.py)
pass builds a local minicram (catalog regions + mates), then EH never touches
GCS again.

This is a **Verily Workbench** adaptation ([`SOURCE.md`](SOURCE.md)). Do **not**
submit it from Terra. In a VWB Jupyter app:

1. Reuse the GATK `Homo_sapiens_assembly38` FASTA + fai already staged by
   [`notebooks/rw/locityper_00_prep_reference.ipynb`](../notebooks/rw/locityper_00_prep_reference.ipynb)
2. [`notebooks/rw/expansion_hunter_01_run.ipynb`](../notebooks/rw/expansion_hunter_01_run.ipynb)
   — stage catalogs + WDL / register / submit via
   [`wb workflow`](https://support.workbench.verily.com/docs/guides/workflows/cromwell/)

The Workbench **Workflows GUI is unreliable**.

## What the WDL does

[`wdl/ExpansionHunterMinicram.wdl`](wdl/ExpansionHunterMinicram.wdl) — workflow
name `ExpansionHunterMinicram`.

| Step | What |
|---|---|
| `MakeMinicram` | Once per sample: `python3 -m str_analysis.make_minicram_for_expansion_hunter` maps the variant catalog through the CRAI, pulls mates, and downloads each CRAM container **once** |
| `ExpansionHunterGenotype` | bw2 `ExpansionHunter --analysis-mode optimized-streaming` on that local minicram |

`minicram_n_preemptible` defaults to **0** so a preempt does not re-read
Nearline. Smoke-test **one sample** and the two-locus catalog before a cohort
CSV.

Do not swap this extract for Locityper’s `print_reads` + BED: EH needs the
mate pass (and any `OfftargetRegions`) so genotypes match a full-CRAM run.

## Run on Verily Workbench

1. Open a JupyterLab app in the workspace (CLI `wb` is on `PATH`).
2. Clone this repo into the app if it is not already a git resource.
3. Build/push the EH image (`./build_docker.sh`) **and** have the print_reads
   image from [`locityper/`](../locityper/) (or mirrors in a workspace-readable
   Artifact Registry). VPC-SC cannot pull the default tags.
4. Confirm `locityper_00_prep_reference` already uploaded assembly38 FASTA +
   fai (EH does not need Jellyfish).
5. Open `expansion_hunter_01_run`. Set `PRINT_READS_DOCKER` / `EH_DOCKER` if
   you mirrored the images. `RECREATE_WORKFLOW=True` once after a WDL change.
6. Stage catalogs + WDL → register `expansion-hunter` → submit. `SUBMIT` is
   off until you set it.

v9 Illumina CRAMs:
`workspace/vwb-aou-datasets-controlled-v9/v9/wgs/cram/manifest.csv`
(`gs://vwb-aou-datasets-controlled/pooled/wgs/cram/v8_base/wgs_{person_id}.cram`).

Staged reference (Locityper prep notebook, 2026-09-15):

```
gs://aou-lr-phase2-resources/locityper/refs/Homo_sapiens_assembly38.fasta
gs://aou-lr-phase2-resources/locityper/refs/Homo_sapiens_assembly38.fasta.fai
```

Smoke catalog (upload from the runner notebook):
[`configs/smoke.catalog.json`](configs/smoke.catalog.json) — first two loci from
the 711-locus panel. Full panel:
[`configs/candidate_EH_Loci.GRCh38.json`](configs/candidate_EH_Loci.GRCh38.json).
This catalog has **no** `OfftargetRegions`; minicram extract is ReferenceRegion
+ mates only.

Single job:

```bash
wb workflow job run \
  --workflow=expansion-hunter \
  --output-bucket-id="$OUTPUT_BUCKET_ID" \
  --output-path=workflowRuns/expansion_hunter \
  --inputs-uri=gs://BUCKET/expansion_hunter/inputs.json
```

Batch (one row per sample). Shared files must be repeated on every row;
Workbench does not fill WDL defaults from a sidecar JSON when you use a CSV:

```bash
wb workflow job run \
  --workflow=expansion-hunter \
  --output-bucket-id="$OUTPUT_BUCKET_ID" \
  --output-path=workflowRuns/expansion_hunter \
  --batch-input-bucket-id="$OUTPUT_BUCKET_ID" \
  --batch-input-csv-path=expansion_hunter/batch.csv \
  --column-mapping-uri=gs://BUCKET/expansion_hunter/column_mapping.json
```

Examples: [`configs/expansion_hunter.inputs.json.example`](configs/expansion_hunter.inputs.json.example),
[`configs/batch.header.csv`](configs/batch.header.csv),
[`configs/column_mapping.json`](configs/column_mapping.json).

Register (after copying the WDL into a workspace bucket):

```bash
wb workflow create \
  --bucket-id="$OUTPUT_BUCKET_ID" \
  --path=wdl/expansion_hunter/ExpansionHunterMinicram.wdl \
  --workflow=expansion-hunter \
  --workflow-type=WDL \
  --display-name="ExpansionHunter minicram"
```

Git repos are **not** valid workflow sources; the WDL has to live in GCS.
See [Create batch jobs](https://support.workbench.verily.com/docs/guides/workflows/create_batch_jobs/).

## Inputs

| Input | Notes |
|---|---|
| `sample_id` | Output prefix |
| `cram` / `crai` | Mapped WGS. CRAM is `localization_optional` (stays in Nearline); CRAI is localized |
| `ref_fa` / `ref_fai` | Uncompressed FASTA + fai (AoU srWGS: GATK `Homo_sapiens_assembly38`) |
| `catalog` | ExpansionHunter variant catalog JSON. Smoke: [`configs/smoke.catalog.json`](configs/smoke.catalog.json) |
| `sex` | `male` / `female` (also `m`/`f`/`1`/`2`). Default `female` |
| `gcloud_project` | Requester-pays billing project for the CRAM bucket; empty uses `$GOOGLE_CLOUD_PROJECT` |
| `print_reads_docker` | Default `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.0` (same image as Locityper `MakeMinicram`) |
| `eh_docker` | Default `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-expansion-hunter:0.1.0` ([`build_docker.sh`](build_docker.sh)) |

VPC-SC workspaces often cannot pull Docker Hub or Broad AR. Mirror
`print_reads_docker` and `eh_docker` to a workspace-readable Artifact Registry.

## Outputs

| Output | What |
|---|---|
| `eh_json` | `{sample_id}.EH.json` |
| `eh_vcf` | `{sample_id}.EH.vcf` |
| `minicram` / `minicrai` | Per-sample subset written by `make_minicram_for_expansion_hunter` |
| `minicram_transfer_stats` | Bytes / containers pulled from Nearline |

Monitor with `wb workflow job list --workflow=expansion-hunter` and
`wb workflow job describe --job-id=...`. Cancel with
`wb workflow job cancel --job-id=...`.
