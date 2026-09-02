# FELIX / Tractor-Mix evaluation gates

Scripts for milestone gates before genome-wide Terra runs. Run from the repo
root with the pilot / genome result TSVs on local disk or synced from Terra.

## Prerequisites

```bash
pip install pandas numpy matplotlib
# optional: scipy (faster erfcinv in calibration scripts)
```

Stage workflow scripts to your workspace bucket first:

```bash
WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_felix_scripts.sh
WORKSPACE_BUCKET=gs://fc-secure-... ./scripts/stage_tractor_scripts.sh
```

## M1 — FELIX vs Tractor-Mix chr22 calibration (2×2)

Compare limited/full covariate runs on the **same cohort** and phenotypes:

| Model            | Covariates | Result TSV suffix   |
|------------------|------------|---------------------|
| `felix_limited`  | limited    | `.felix.tsv`        |
| `felix_full`     | full       | `.felix.tsv`        |
| `tractor_limited`| limited    | `.tractor_mix.tsv`  |
| `tractor_full`   | full       | `.tractor_mix.tsv`  |

```bash
python3 felix/scripts/compare_felix_tractor_calibration.py \
  --felix-limited-dir /path/to/felix_chr22_limited/results_by_phenotype \
  --felix-full-dir    /path/to/felix_chr22_full/results_by_phenotype \
  --tractor-limited-dir /path/to/tractor_chr22_limited/results_by_phenotype \
  --tractor-full-dir    /path/to/tractor_chr22_full/results_by_phenotype \
  --out-dir eval/chr22_felix_tractor_gate
```

**Pass criteria:** λGC ≈ 1 on null phenotypes for all four models; `comparability_flags.tsv`
status `OK` (matched effective N). Do not interpret λGC differences when status is `FAIL`.

Outputs:

- `calibration_summary.tsv` — per-phenotype λGC and tested-variant counts
- `lambda_gc_wide.tsv` — wide pivot for notebooks
- `<phenotype>/qq_matched.png` — overlay QQ for the four models

## M4 — Repeat-unit encoding comparison

After `FelixGenome.wdl` RU_TEST branch (or local `tractor-mix-score --mode felix` runs):

```bash
python3 felix/scripts/compare_repeat_encodings.py \
  --results-dir /path/to/ru_scores/chr22 \
  --phenotype my_trait \
  --out-dir eval/chr22_encoding_compare
```

Expects files named like `my_trait.chr22.ru_dosage.felix.tsv`,
`...ru_collapse.felix.tsv`, `...ru_split.felix.tsv`.

Outputs: `calibration_summary.tsv`, `encoding_comparison.tsv`, per-encoding QQ plots.

## M4 — Repeat dosage simulator

Quick type I error / power check for encoding architectures (no VCF required):

```bash
python3 felix/scripts/simulate_repeat_dosage.py \
  --n-samples 2000 \
  --n-reps 500 \
  --effect 0.15 \
  --out-dir eval/sim_repeat_dosage
```

Architectures: `length_additive`, `single_allele`, `threshold`.

Outputs:

- `simulation_summary.tsv` — reject rates at `--alpha` (default 0.05)
- `type1_error.tsv` — null-effect rows only
- `power.tsv` — alternative-effect rows only

## Terra workflows

| WDL | Scope |
|-----|-------|
| `felix/wdl/FelixPilot.wdl` | chr22 pilot (felixla + step2) |
| `felix/wdl/FelixGenome.wdl` | Autosomes: shared null, per-chr pack/step2, optional RU_TEST scorer |

Example inputs: `felix/configs/felix.genome.inputs.{limited,full}.json.example`

Docker default: `us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0`

Build / push image:

```bash
felix/build_docker.sh
```

## Notes

- FELIX joint p-value for calibration: `P_cct_admixed_c`
- Tractor-Mix joint p-value: `P`
- Genome-wide hit tables come from `felix/scripts/summarize_felix_results.py` inside the WDL summarize task
- Do not headline encoding “more hits” without matched λGC
