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

1. Build / push the Docker image if needed (tag `0.3.2`):

```bash
cd tractor_mix
./build_docker.sh
```

2. Stage scripts + rebuilt covariates on the workspace bucket:

```bash
gsutil -m cp scripts/* "$WORKSPACE_BUCKET/scripts/"
gsutil cp tractor_mix/covariates.source_rebuilt.csv.gz \
  "$WORKSPACE_BUCKET/covariates/covariates.source_rebuilt.csv.gz"
```

3. Run `notebooks/tractor_01_prepare_inputs.ipynb` on the Workbench.
   It pulls `covariates.source_rebuilt.csv.gz` (default URI above, or
   `TRACTOR_COVARIATES_GCS`), builds the shared recommended-model-complete
   cohort, and uploads under `$WORKSPACE_BUCKET/tractor_mix_pilot/`:
   `analysis_samples.txt`, `pheno_cov.tsv`, `covariate_columns_{limited,full}.txt`,
   etc. Confirm `analysis_samples` line count ≫ 1000 before submitting.

4. Submit **Tractor-Mix only** first (same cohort; only `covariate_columns` differs):

   - `TractorMixPilot.wdl` + `configs/pilot.inputs.limited.json.example`
   - `TractorMixPilot.wdl` + `configs/pilot.inputs.full.json.example`

   Point `gs://BUCKET/...` at your workspace bucket. WDL script inputs should
   use `$WORKSPACE_BUCKET/scripts/`.

5. After Tractor-Mix succeeds, optionally submit SAIGE with the same cohort
   (`configs/saige.inputs.*.json.example`), then run
   `notebooks/tractor_02_qc_results.ipynb`.

## Cohort covariates (rebuild + atlas)

1. Place CDRv9 source exports in `resources/` (gitignored).
2. Run `notebooks/tractor_00_cov_rebuild_source.ipynb` →
   `covariates.source_rebuilt.csv.gz` (+ data dictionary).
3. Run `notebooks/tractor_03_cov_summarize.ipynb` → inline QC plus
   `summaries/{figures,tables,manuscript}/` (gitignored). Small crosstab cells
   (`n < 20`) are redacted in exports.

## Notes

- Tractor-Mix: unconditional `TractorMix.score`, `AC_threshold=50`, logistic nulls.
- FLARE extract uses the local Rust `extract-tracts-flare` binary (not
  `/opt/Tractor/scripts/extract_tracts_flare.py`). See
  `extract_tracts_flare/README.md`.
- SAIGE: sparse GRM from the same chr1+chr22 markers, logistic null + SPA step2.
- Global PCs from covariates (not PC-AiR); relatedness cutoff default `0.05`.
- Image `0.3.2` ships an Alpine-static `extract-tracts-flare` and checks that
  it can exec on `wzhou88/saige:1.3.3` at image-build time.
