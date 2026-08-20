# tractor-mix-score — upstream parity and performance notes

This document is written for [Atkinson-Lab Tractor-Mix](https://github.com/Atkinson-Lab/Tractor-Mix)
maintainers and anyone evaluating a Rust reimplementation of the scoring step.

**Claim:** `tractor-mix-score` is a drop-in replacement for the **sparse GRM** code path in
`TractorMix.score.R`. It does **not** change association statistics at the precision Tractor-Mix
already reports. The goal is a large speedup on whole-genome scoring, not a new method.

## Reference implementation

| Item | Value |
|------|-------|
| Oracle script | [`oracle/TractorMix.score.R`](oracle/TractorMix.score.R) |
| Upstream commit | [`4adb8f18`](https://github.com/Atkinson-Lab/Tractor-Mix/blob/4adb8f1814d9315ecd7868eb729d52ec0c723719/TractorMix.score.R) |
| Branch / path | Sparse GRM (`Sigma_i`, `Sigma_iX`) only |
| Null model | Unchanged — still fit in R with `GMMAT::glmmkin` |

We vendor the R script at a fixed SHA so parity tests always compare against a known baseline.

## Scope (intentional limits)

- **Binomial** nulls only (`family = binomial(link = "logit")`).
- **Sparse kinship** path only (exported `Sigma_i` as CSC matrix).
- Does **not** reimplement `glmmkin`, ridge shrinkage, or null-model fitting.
- Does **not** support dense-only or non-sparse code paths in upstream `TractorMix.score.R`.

These match our All of Us long-read pilot requirements. Extending scope for upstream would be a
separate discussion.

## Architecture

```
R (unchanged)                    Rust (new)
─────────────                    ──────────
glmmkin → null.rds               —
export CSC binaries ────────────→ load null_export/
TractorMix.score.R (oracle)      tractor-mix-score (production)
     ↓                                    ↓
  phenotype.tsv                    phenotype.tsv  (same schema)
```

R exports a small binary bundle after null fit (`scripts/fit_null.R`):

- `sigma_i.csc.bin` — sparse `(Sigma_i)` in CSC layout
- `sigma_i_x.bin`, `cov.bin`, `residuals.bin`, `id_include.txt`, `meta.json`

Rust streams gzip dosage files and writes the same TSV columns as Tractor-Mix.

## Parity methodology

1. **Synthetic null** — `glmmkin` on a random sparse kinship matrix (`testdata/gen_oracle.R` or
   `testdata/bench_large_fixture.R`).
2. **Synthetic dosages** — multi-ancestry `.dosage.txt.gz` files with controlled AC patterns
   (all pass filter, partial drop, all drop, joint NA cases).
3. **R oracle** — vendored `TractorMix.score.R` via `testdata/run_oracle.R`.
4. **Rust scorer** — `tractor-mix-score` on the exported null bundle.
5. **Compare** — `testdata/compare_oracle_tsv.py` applies the same post-processing Tractor-Mix uses:
   - `round(x, 5)` for `Chi2`, `Eff_anc*`, `SE_anc*`
   - `signif(x, 5)` for `P`, `Pval_anc*`
   - exact match for variant IDs and integer columns

### Reproduce parity (recommended)

Inside the pilot Docker image (GMMAT + prebuilt binaries):

```bash
cd tractor_mix/tractor_mix_score
./scripts/run_oracle_parity_in_docker.sh
```

Expected output:

```
OK: 48 rows match between .../r.tsv and .../rust.tsv
Oracle parity passed.
```

On this **official 48-variant fixture**, R and Rust outputs are **byte-identical** TSV files —
not merely “close” at floating-point tolerance.

### Pre-Terra performance gate (required)

```bash
./scripts/bench_realistic_in_docker.sh
```

Compares **R `TractorMix.score` with `n_core=4`** (legacy production default) vs **Rust
`--threads 8`** on a thousands-of-samples × tens-of-thousands-of-variants fixture.
**Do not submit Terra scoring until this passes** (Rust faster than parallel R).

Toy benchmark (`./scripts/bench_r_vs_rust_in_docker.sh`, 800 × 2000 variants) is for smoke
tests only — do not extrapolate to chr22 / ~10k samples.

### Toy benchmark (smoke test only)

| Step | Wall time (example) |
|------|---------------------|
| R `TractorMix.score.R` (4 cores) | varies |
| Rust `tractor-mix-score` (8 threads) | varies |

At AoU scale, null fit (`glmmkin`) remains ~3h/phenotype in R; only the **Score** task shrinks.

### Terra chr22 limited pilot (Aug 2026)

Real cohort run with `tractor-mix-pilot:0.4.2`, `--threads 8`:

| Metric | Value |
|--------|-------|
| Samples (`n`) | 9,207 |
| Variants (chr22) | 169,213 |
| Ancestries | 5 |
| Score throughput | ~98 variants/sec |
| Score wall time | ~29 min compute (~34–45 min shard including overhead) |

Compared with the prior R combined FitNull+Score (~7h/phenotype) and single-thread Rust 0.4.1 Score (~1.5h+), parallel Rust Score is the production path.

## What changed vs what did not

| Unchanged | Changed |
|-----------|---------|
| Null model fitting (`glmmkin`) | Scoring engine implementation language |
| Input dosage file format | Null objects passed as exported CSC binaries |
| Output TSV schema and R formatting rules | Score wall time (when `--threads` used and realistic bench passes) |
| AC filter logic, joint test, per-ancestry tests | Parallel variant scoring within chunks (`rayon`, `--threads`) |

## Suggested upstream contribution shape

If Atkinson-Lab is interested, a minimal upstream-friendly package could include:

1. **`tractor-mix-score/`** — Rust crate + CLI (this directory).
2. **`R/export_null_for_score.R`** — companion to emit the CSC bundle from an existing
   `glmmkin` object (our `scripts/fit_null.R` export block).
3. **Parity tests** — Docker script + vendored oracle; no GMMAT required for Rust unit tests.
4. **This document** — parity evidence and benchmark instructions.

We are happy to adapt naming, install story (CRAN vs standalone binary), and scope to match
upstream preferences.

## Local development (no Docker)

```bash
cd tractor_mix/tractor_mix_score
cargo test                                    # R-free unit tests
cargo test -- --ignored --nocapture           # full oracle (needs local GMMAT)
```

## Contact

Developed for the All of Us long-read phase-2 Tractor-Mix pilot (Broad LRMA). For questions about
parity fixtures, benchmarks, or contributing back to Tractor-Mix, open an issue or PR in the
repository that contains this crate.
