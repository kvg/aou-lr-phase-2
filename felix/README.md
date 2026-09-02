# FELIX LAI GWAS (Terra)

Chr22 pilot and genome-wide workflows for FELIX local-ancestry-aware association on
All of Us / Terra. Consumes the published `lhu1/felix:latest` image (do not compile
or fork FELIX). Shares cohort inputs with Tractor-Mix / SAIGE pilots under
`tractor_mix_pilot/`.

## Layout

```
felix/
  build_docker.sh              # Cloud Build or local docker (repo-root context)
  docker/Dockerfile            # FROM lhu1/felix:latest + PLINK2/bcftools + Rust CLIs
  docker/cloudbuild.yaml
  wdl/FelixPilot.wdl           # chr22: mtx GRM + felixla + admixed step2
  wdl/FelixGenome.wdl          # autosomes: shared null, per-chr pack/step2, optional RU_TEST scorer
  configs/felix*.json.example
  scripts/                     # FELIX-specific CLIs (staged to gs://BUCKET/felix/scripts/)
  eval/README.md               # M1/M4 calibration gates
scripts/                       # shared: GRM build, make_plink_keep, annotate_repeat_units, …
tractor_mix/
  extract_tracts_flare/        # RU dosage extraction (used by FelixGenome RU_TEST branch)
  tractor_mix_score/           # Rust scorer (--mode felix)
```

## Design

FELIX replaces Tractor-Mix for SNV/indel + unique-SV association (`felixla` pack +
`step2_SPAtests` with `--is_admixed=TRUE`). Repeat-mediated SVs use the optional
**RU_TEST** branch: annotate → `extract-tracts-flare` (dosage encoding) →
`tractor-mix-score --mode felix` with the same FELIX null export.

| Model | Covariates (2×2 calibration) |
|-------|------------------------------|
| FELIX limited | sex, age, PC1–PC10 |
| FELIX full | limited + standardized coverage + GC dummies |

Scan: chr22 FLARE / phased VCFs (pilot) or parallel `chroms` + `phase_vcfs` +
`flare_vcfs` (genome). GRM markers: chr1 + chr22 (same as Tractor-Mix / SAIGE).

## Build image

```bash
felix/build_docker.sh                 # Cloud Build → us-central1-docker.pkg.dev/.../felix-pilot:0.1.0
felix/build_docker.sh --local --push  # or local docker
```

Bump `FelixPilot.docker` / `FelixGenome.docker` in config JSON when retagging.

## Stage scripts

FELIX scripts live under `felix/scripts/`; shared helpers stay in repo `scripts/`:

```bash
WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_felix_scripts.sh
WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_tractor_scripts.sh   # GRM, annotate_repeat_units, …
```

## Run order

1. Prepare cohort with `notebooks/terra/tractor_01_prepare_inputs.ipynb` (same as Tractor-Mix).
2. Submit **FelixPilot.wdl** with `felix/configs/felix.pilot.inputs.{limited,full}.json.example`.
   Point `gs://BUCKET/...` at your workspace bucket.

   Workflow: **Check** → **MakeGRM** (`build_saige_plink_and_grm.sh`) → **FitFelixNull**
   (`fit_felix_null.R` → FELIX step1 + `export_felix_null.R`) → **Pack** (`felixla`) →
   **Step2** (`run_felix_step2.R` / `step2_SPAtests.R`, `--is_admixed=TRUE`) → summarize
   (`summarize_felix_results.py`, p-column `P_cct_admixed_c`).

3. Compare λGC vs Tractor-Mix 2×2: `felix/eval/README.md`.

4. Genome-wide: **FelixGenome.wdl** + `felix/configs/felix.genome.inputs.*.json.example`.
   Optional RU_TEST: set `simple_repeat_bed` + `annotate_repeat_units_script`
   (`gs://BUCKET/scripts/annotate_repeat_units.py`).

Resolve FLARE URIs:

```bash
python3 scripts/resolve_flare_uris.py --from-firecloud --scan-chrom chr22 --grm-chroms chr1 chr22
python3 scripts/resolve_flare_uris.py --from-firecloud --autosomes --grm-chroms chr1 chr22
```

## Notes

- Joint p-value for calibration: `P_cct_admixed_c` (FELIX) vs `P` (Tractor-Mix).
- `fit_null.R --step1-rda` delegates to `felix/scripts/export_felix_null.R` for local dev.
- Do not headline encoding comparisons without matched λGC (`felix/eval/README.md`).
