# Tractor-Mix + SAIGE Terra calibration pilots

End-to-end chr22 pilots for matched local-ancestry-aware (Tractor-Mix) and
standard mixed-model (SAIGE) association on All of Us / Terra.

## Layout

```
notebooks/                       # all Workbench notebooks (upload with scripts/)
scripts/                         # all CLIs used by notebooks + WDLs
tractor_mix/
  docker/Dockerfile
  extract_tracts_flare/
  tractor_mix_score/
  wdl/TractorMixPilot.wdl
  wdl/SaigePilot.wdl
  configs/pilot.inputs.{limited,full}.json.example
  configs/saige.inputs.{limited,full}.json.example
  resources/                     # local source exports (gitignored)
  summaries/                     # figures + TSV tables from tractor_03_ (gitignored)
sv_annotation/
  wdl/AnnotateSvCallset.wdl
  docker/
  configs/
  tests/
```

## Design (2×2 calibration)

One shared **full-covariate-complete** cohort and phenotype table:

| Model | Covariates |
|-------|------------|
| Tractor-Mix limited | sex, age, PC1–PC10 |
| Tractor-Mix full | limited + standardized coverage + GC dummies |
| SAIGE limited | same limited list |
| SAIGE full | same full list |

Scan: chr22 FLARE VCF. GRM markers: chr1 + chr22.

## FLARE VCF URIs (chr1 + chr22)

Resolved live from the Terra data table `aou_lr_chrom` in workspace
`allofus-drc-wgs-LR-prodData` / `AoU_DRC_LongReads_PhaseTwo_Storage` via firecloud
(`model_chr_anc_vcf`). Column `global_anc` is `*.global.anc.gz`, not the LAI VCF.

```bash
python3 ../scripts/resolve_flare_uris.py --from-firecloud \
  --scan-chrom chr22 --grm-chroms chr1 chr22
```

FLARE ancestries: `eas=0,amr=1,eur=2,afr=3,sas=4` → `num_ancs=5`.

## Run order (Tractor-Mix first)

1. Build / push the Docker image if needed (tag `0.4.2`):

```bash
cd tractor_mix
./build_docker.sh
```

Before Terra scoring, run the realistic benchmark locally or in Docker:

```bash
cd tractor_mix/tractor_mix_score
./scripts/bench_realistic_in_docker.sh   # requires 0.4.2 image
```

2. Stage scripts + rebuilt covariates on the workspace bucket:

```bash
# from repo root — uploads fit_null.R and other WDL CLIs
WORKSPACE_BUCKET=gs://... ./scripts/stage_tractor_scripts.sh

gsutil cp tractor_mix/covariates.source_rebuilt.csv.gz \
  "$WORKSPACE_BUCKET/covariates/covariates.source_rebuilt.csv.gz"
```

Or sync the full scripts tree: `gsutil -m rsync -r scripts/ "$WORKSPACE_BUCKET/scripts/"`

3. Run `notebooks/tractor_01_prepare_inputs.ipynb` on the Workbench.
   It pulls `covariates.source_rebuilt.csv.gz`, builds the shared cohort, uploads
   under `$WORKSPACE_BUCKET/tractor_mix_pilot/`, and can re-stage WDL scripts
   (`fit_null.R`, etc.). Confirm `analysis_samples` line count ≫ 1000 before submitting.

4. Submit **Tractor-Mix only** first (same cohort; only `covariate_columns` differs):

   - `TractorMixPilot.wdl` + `configs/pilot.inputs.limited.json.example`
   - `TractorMixPilot.wdl` + `configs/pilot.inputs.full.json.example`

   Point `gs://BUCKET/...` in the input JSON at your workspace bucket. Use
   `$WORKSPACE_BUCKET/scripts/fit_null.R` (not legacy `fit_null_and_score.R`).

   Workflow: **FitNull** (`fit_null.R` + `glmmkin`, ~3h/phenotype) → **Score**
   (`tractor-mix-score --threads 8`, ~30–45 min/phenotype on chr22). Shared
   `docker` default is `0.4.2`. For Score-only image iteration, override nested
   `TractorMixPilot.Score.docker` (`allowNestedInputs: true`).

5. After Tractor-Mix succeeds, submit SAIGE with the same cohort
   (`SaigePilot.wdl` + `configs/saige.inputs.*.json.example`), then run
   `notebooks/tractor_02_qc_results.ipynb`.

   SAIGE inputs mirror Tractor: same `analysis_samples` / `pheno_cov` /
   `selected_phenotypes` / `grm_vcfs`, plus staged
   `build_saige_plink_and_grm.sh`, `make_plink_keep.py`, `fit_saige_null.R`,
   `run_saige_step2.R`. Default `min_mac=50` matches Tractor `ac_threshold`.
   Image `0.4.2` already includes SAIGE 1.3.3, PLINK2, and bcftools — no
   rebuild needed for SAIGE-only iteration.

## Cohort covariates (rebuild + atlas)

1. Place CDRv9 source exports in `resources/` (gitignored).
2. Run `notebooks/tractor_00_cov_rebuild_source.ipynb` →
   `covariates.source_rebuilt.csv.gz` (+ data dictionary).
3. Run `notebooks/tractor_03_cov_summarize.ipynb` → inline QC plus
   `summaries/{figures,tables,manuscript}/` (gitignored). Small crosstab cells
   (`n < 20`) are redacted in exports.

## Notes

- Tractor-Mix: Rust `tractor-mix-score` (sparse GRM branch of pinned `TractorMix.score`), `AC_threshold=50`, logistic nulls via `fit_null.R` + `glmmkin`.
- FLARE extract uses the local Rust `extract-tracts-flare` binary (not
  `/opt/Tractor/scripts/extract_tracts_flare.py`). See
  `extract_tracts_flare/README.md`.
- Score tests use `tractor_mix_score/` (`tractor-mix-score` CLI); null fit stays R.
  R-oracle parity and benchmarks: `./scripts/run_oracle_parity_in_docker.sh`,
  `./scripts/bench_realistic_in_docker.sh`. Upstream contribution notes:
  `tractor_mix_score/UPSTREAM.md`.
- SAIGE: sparse GRM from the same chr1+chr22 markers, logistic null + SPA step2
  (`--vcfField=GT`, `min_mac=50`). BuildGRM / FitNull use `preemptible=0`.
- Global PCs from covariates (not PC-AiR); relatedness cutoff default `0.05`.
- Image `0.4.2` ships Alpine-static `extract-tracts-flare` and `tractor-mix-score` (parallel `--threads`), verified to exec on `wzhou88/saige:1.3.3` at image-build time.
