# Terra notebooks

Reviewer-facing map (manuscript tables, run order, env vars):
[`../../README.md`](../../README.md).

Terra analysis notebooks live here. Companion CLIs live in [`../../scripts/`](../../scripts/).

## Start here on Terra: sync from GitHub

AoU Research Program workspaces are not reachable with laptop `gsutil`. Run
everything from a Terra Workbench notebook (pmi-ops browser login).

1. **Once:** upload [`00_sync_repo.ipynb`](00_sync_repo.ipynb) into the workspace
   (or paste its cells).
2. Run all cells (repo is public HTTPS — no token). That will:
   - `git clone` / `git pull` the repo on the VM
   - `gsutil rsync` `scripts/` → `$WORKSPACE_BUCKET/scripts/`
   - `gsutil rsync` terra notebooks → `$WORKSPACE_BUCKET/notebooks/`
   - stage WDLs under `$WORKSPACE_BUCKET/wdl/`
   - push a new `FlareByPopulation` method snapshot + bump the workspace config
   - upsert `flare_lai_exp` via FISS
   - optionally submit incomplete `flare_lai_exp` rows (`SUBMIT_FLARE=True`)
3. Re-run `00_sync_repo` whenever you need the latest `main` on Terra.

Helper CLI (same VM): [`../../scripts/terra_sync_repo.py`](../../scripts/terra_sync_repo.py).

## Script bootstrap on Terra

On a Terra notebook VM, pipeline CLIs live at `$WORKSPACE_BUCKET/scripts/`. Persistent `edit/scripts/` copies are often stale and must not win. Each notebook’s first code cell:

1. `gsutil -m rsync` `$WORKSPACE_BUCKET/scripts/` into `./scripts/` when `WORKSPACE_BUCKET` is set (`TERRA_SYNC_SCRIPTS=false` skips this).
2. Drops cached `sys.modules` entries from that directory so a kernel re-run sees the new files.
3. Calls `terra_notebook.init_notebook(...)` and verifies required CLIs.

Prefer **`00_sync_repo`** over laptop staging. If you already have a clone on the
VM, you can still run:

```bash
python3 aou-lr-phase-2/scripts/terra_sync_repo.py --ref main
```

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

SNV / indel callset stats:

| Notebook | Use |
|---|---|
| `snv_00_merge_glnexus_stats.ipynb` | Fetch `GL_INTERVAL_set.stats`, merge Table 2 SNV/indel counts |

Methylation maps:

| Notebook | Use |
|---|---|
| `meth_00_merge_pbcpg_stats.ipynb` | Firecloud-fetch `aou2_v1_phased_bams` from the storage workspace, pull `stats_tsv`, merge, optional concordance |

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
| `00_sync_repo.ipynb` | **Start here:** clone/pull repo on the Terra VM; stage scripts, notebooks, WDLs; upsert tables via FISS |
| `flare_01_switch_gq_dp.ipynb` | Ancestry-switch GQ/DP QC, plus old vs new FLARE switch rates / tract lengths |
| `flare_02_lai_exp_compare.ipynb` | Firecloud-fetch `flare_lai_exp` method grid; switch-QC finished rows as they complete |
| `flare_03_stage_site_filters.ipynb` | Stage Part 4 context mask + Part 7 call-QC `include_sites` (+ optional concordance) |

Per-population FLARE (own T per `population`): `flare/wdl/FlareByPopulation.wdl`. See `flare/README.md`.
Copy covering FLARE `AN1`/`AN2` onto interstitial sites of a phased target VCF:
`propagate_annotations/wdl/PropagateFlareAncestry.wdl`.

WDLs stay in `tractor_mix/wdl/`, `felix/wdl/`, `sv_annotation/wdl/`,
`propagate_annotations/wdl/`, `snv_stats/wdl/`, `methylation_stats/wdl/`,
and `flare/wdl/`.

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
