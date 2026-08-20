# Workbench notebooks

All analysis notebooks live here. Companion CLIs live in [`../scripts/`](../scripts/).

Upload **both** `notebooks/` and `scripts/` into the Terra workspace (siblings). Jupyter’s working directory can be either the workspace root or `notebooks/`.

## Script bootstrap on Terra

On a Terra notebook VM, pipeline CLIs usually live only at `$WORKSPACE_BUCKET/scripts/`. Each notebook’s first code cell calls `terra_notebook.init_notebook(...)` to:

1. Use a local sibling `scripts/` directory when present (git checkout or prior sync).
2. Otherwise `gsutil rsync` the full `$WORKSPACE_BUCKET/scripts/` prefix into `./scripts/`.
3. Verify required files exist (fetching individual missing objects if needed).

Upload `scripts/terra_notebook.py` and `scripts/workspace_paths.py` with the rest of `scripts/` before running notebooks on Terra.

Tractor-Mix (run in this order):

| Notebook | Use |
|---|---|
| `tractor_00_cov_rebuild_source.ipynb` | Rebuild person-level covariates |
| `tractor_01_prepare_inputs.ipynb` | Tractor-Mix / SAIGE pilot inputs |
| `tractor_02_qc_results.ipynb` | Calibration QC |
| `tractor_03_cov_summarize.ipynb` | Covariate atlas |
| `tractor_04_table1_cohort_summary.ipynb` | Manuscript Table 1 |
| `tractor_05_pca_deepvariant_long_read.ipynb` | Hail DeepVariant PCA |
| `tractor_06_qc_grm.ipynb` | GRM kinship QC |

SV annotation (run in this order):

| Notebook | Use |
|---|---|
| `sv_00_prepare_caddsv_annotations.ipynb` | One-time CADD-SV bundle (≥60 GB disk) |
| `sv_01_stage_repeat_tracks.ipynb` | UCSC rmsk / simpleRepeat / genomicSuperDups BEDs |
| `sv_02_stage_sample_ancestry.ipynb` | Phase 1 / Phase 2 sample ancestry TSVs |
| `sv_03_manuscript_stats.ipynb` | SV site-count table + discovery plots (after the WDL) |

WDLs stay in `tractor_mix/wdl/` and `sv_annotation/wdl/`.

## Tractor-Mix WDL (Rust scorer)

`TractorMixPilot.wdl` runs **FitNull** (`scripts/fit_null.R`) then **Score** (`tractor-mix-score --threads 8`).
Shared `docker` default is `tractor-mix-pilot:0.4.2`. For Score-only image iteration, override nested
`TractorMixPilot.Score.docker` (`allowNestedInputs: true`).
Run `tractor_mix/tractor_mix_score/scripts/bench_realistic_in_docker.sh` before production Terra scoring.

Before submitting on Terra, stage scripts (must include `fit_null.R`):

```bash
WORKSPACE_BUCKET=gs://... ./scripts/stage_tractor_scripts.sh
```

Do **not** use legacy `fit_null_and_score.R` on the workspace bucket.
