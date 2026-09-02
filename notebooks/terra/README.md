# Terra notebooks

Terra analysis notebooks live here. Companion CLIs live in [`../../scripts/`](../../scripts/).

Upload the contents of `notebooks/terra/` and `scripts/` into the Terra workspace as siblings. Jupyter’s working directory can be either the workspace root or the notebooks folder.

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
| `tractor_05_pca_deepvariant_long_read.ipynb` | Hail DeepVariant **global** PCA (genotype-only; no covariates) |
| `tractor_06_fill_lr_ancestry.ipynb` | Fill remaining long-read ancestry gaps (`lr_PC` kNN) |
| `tractor_07_pca_within_population.ipynb` | Hail within-population PCA (after filled covariates + `qc_for_pca.mt`) |
| `tractor_08_qc_grm.ipynb` | GRM kinship QC |
| `tractor_09_genome_post_workflow.ipynb` | TractorMixGenome post-workflow QC / PheWAS views |

Both PCA notebooks (`tractor_05` and `tractor_07`) are **resumable**: completed stages are skipped when
their outputs already exist under the same `PCA_RUN_LABEL`. Use `PCA_FORCE_REIMPORT`,
`PCA_FORCE_GLOBAL_PCA`, or `PCA_FORCE_WITHIN_POPULATION` to rerun a stage. See each
notebook’s intro for details.

SV annotation (run in this order):

| Notebook | Use |
|---|---|
| `sv_00_prepare_caddsv_annotations.ipynb` | One-time CADD-SV bundle (≥60 GB disk) |
| `sv_01_stage_repeat_tracks.ipynb` | UCSC rmsk / simpleRepeat / genomicSuperDups BEDs |
| `sv_02_stage_sample_ancestry.ipynb` | Phase 1 / Phase 2 sample ancestry TSVs |
| `sv_03_manuscript_stats.ipynb` | SV site-count table + discovery plots (after the WDL) |

FLARE annotation QC:

| Notebook | Use |
|---|---|
| `flare_01_switch_gq_dp.ipynb` | Ancestry-switch sites vs `GQ`/`DP` (after `PropagateAnnotations.wdl`) |

WDLs stay in `tractor_mix/wdl/`, `felix/wdl/`, `sv_annotation/wdl/`, and `propagate_annotations/wdl/`.

FELIX workflows (`felix/wdl/FelixPilot.wdl`, `FelixGenome.wdl`) use the
`felix-pilot:0.1.0` image. Stage FELIX scripts with `./scripts/stage_felix_scripts.sh`
and shared helpers with `./scripts/stage_tractor_scripts.sh`. See `felix/README.md`.

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
