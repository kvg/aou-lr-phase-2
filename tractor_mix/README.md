# Tractor-Mix + SAIGE Terra calibration pilots

End-to-end chr22 pilots for matched local-ancestry-aware (Tractor-Mix) and
standard mixed-model (SAIGE) association on All of Us / Terra.

## Layout

```
notebooks/terra/                 # Terra notebooks (upload with scripts/)
notebooks/rw/                    # Verily Workbench notebooks
scripts/                         # all CLIs used by notebooks + WDLs
tractor_mix/
  docker/Dockerfile
  extract_tracts_flare/
  tractor_mix_score/
  wdl/TractorMixPilot.wdl        # chr22 calibration pilot
  wdl/TractorMixGenome.wdl       # multi-chr: shared null, scatter extract/score
  wdl/SaigePilot.wdl
  configs/pilot.inputs.{limited,full}.json.example
  configs/genome.inputs.limited.json.example
  configs/saige.inputs.{limited,full}.json.example
  resources/                     # local source exports (gitignored)
  summaries/                     # figures + TSV tables from tractor_03_ (gitignored)
felix/                           # FELIX LAI GWAS (see felix/README.md)
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
# Pilot (single scan chrom)
python3 ../scripts/resolve_flare_uris.py --from-firecloud \
  --scan-chrom chr22 --grm-chroms chr1 chr22

# Genome-wide TractorMixGenome inputs (chroms + flare_vcfs + grm_vcfs)
python3 ../scripts/resolve_flare_uris.py --from-firecloud \
  --autosomes --grm-chroms chr1 chr22
```

FLARE ancestries: `eas=0,amr=1,eur=2,afr=3,sas=4` → `num_ancs=5`.

### Genome-wide (`TractorMixGenome.wdl`)

Same cohort / scripts / docker as the pilot. Differences:

- `flare_vcfs` + parallel `chroms` (Terra chrom-set friendly)
- **MakeGRM** + **FitNull** once; **Extract** / **Score** scatter over chromosomes
- **Concat** merges chrom shards → one `{phenotype}.tractor_mix.tsv` per phenotype
  (`results_tsvs`); per-chrom shards kept as `results_tsvs_by_chrom`
- **Summarize** copies phenotype-named TSVs (`results_tsvs_named`), writes
  `results_manifest.tsv`, λGC / QQ / Manhattan / top hits, and
  `phewas_genomewide_hits.tsv`
- `grm_vcfs` stays separate and small by default (chr1+chr22) — do not feed all autosomes into GRM on the first pass
- Larger Extract/MakeGRM disks scale with `size(vcf)` (plus a floor); Check only
  compares array lengths and does **not** localize FLARE VCFs.

Post-workflow notebook: `notebooks/terra/tractor_09_genome_post_workflow.ipynb`
(cross-phenotype λ / hit plots, optional limited-vs-full comparison).

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

3. Run `notebooks/terra/tractor_01_prepare_inputs.ipynb` on the Workbench.
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
   `notebooks/terra/tractor_02_qc_results.ipynb`.

   SAIGE inputs mirror Tractor: same `analysis_samples` / `pheno_cov` /
   `selected_phenotypes` / `grm_vcfs`, plus staged
   `build_saige_plink_and_grm.sh`, `make_plink_keep.py`, `fit_saige_null.R`,
   `run_saige_step2.R`. Default `min_mac=50` matches Tractor `ac_threshold`.
   Image `0.4.2` already includes SAIGE 1.3.3, PLINK2, and bcftools — no
   rebuild needed for SAIGE-only iteration.

## Cohort covariates (rebuild + atlas)

1. Place CDRv9 source exports in `resources/` (gitignored).
2. Run `notebooks/terra/tractor_00_cov_rebuild_source.ipynb` →
   `covariates.source_rebuilt.csv.gz` (+ data dictionary).
3. Merge long-read global PCs + HPRC/HGSVC3 control rows (re-run after each
   rebuild):

```bash
python3 scripts/merge_lr_global_pcs_into_covariates.py \
  --global-pcs tractor_mix/pca/deepvariant_lr_v1/global_pcs.tsv \
  --control-metadata tractor_mix/reference_controls/control_sample_metadata.tsv
```

   This adds `lr_PC1`–`lr_PC32`, `has_lr_pcs`, `is_reference_control`, and 293
   `HG*`/`NA*` control rows (ancestry/sex/SV fills where known).
4. Fill remaining long-read ancestry gaps via population copy + `lr_PC` kNN:
   `notebooks/terra/tractor_06_fill_lr_ancestry.ipynb` (writes audit TSV under
   `reference_controls/`). Required before within-pop PCA.
5. After `tractor_07_pca_within_population.ipynb`, merge within-pop PCs:

```bash
python3 scripts/merge_lr_pop_pcs_into_covariates.py \
  --population-pcs tractor_mix/pca/deepvariant_lr_v1/population_pcs.tsv
```

   Adds `lr_pop_PC1`–`lr_pop_PC32`, `has_lr_pop_pcs`, and `lr_pop_population`
   (does not replace short-read `pop_PC*`).
6. Soft joint-callset nits (`population` / `sex_at_birth`):

```bash
python3 scripts/apply_lr_soft_field_fills.py \
  --audit tractor_mix/reference_controls/lr_soft_field_fills.tsv
```

   Replays the audit TSV. Add `--discover` after a rebuild to catch any new
   `has_lr_pcs` gaps with the same rules (`lr_pop_population` → `population`,
   `inferred_sex` XX/XY → Female/Male). Does not overwrite sentinels like
   `PMI: Skip`.
7. Run `notebooks/terra/tractor_03_cov_summarize.ipynb` → inline QC plus
   `summaries/{figures,tables,manuscript}/` (gitignored). Small crosstab cells
   (`n < 20`) are redacted in exports.

Fill provenance (controls, PCs, kNN ancestry, soft fields):
`reference_controls/COVARIATE_FILLS.md`.

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
