# tractor-mix-score

Streaming Rust replacement for Atkinson-Lab Tractor-Mix
[`TractorMix.score.R`](https://github.com/Atkinson-Lab/Tractor-Mix/blob/4adb8f1814d9315ecd7868eb729d52ec0c723719/TractorMix.score.R)
(**sparse GRM branch only**, SHA `4adb8f18`). Binomial nulls only; `glmmkin` stays in R.

**Same results as R; faster scoring at production scale.** Output matches the pinned R
oracle on the parity fixture (byte-identical TSV). `--threads` parallelizes variant scoring
within each chunk (deterministic output order).

**Do not extrapolate from toy benchmarks.** Validate with `./scripts/bench_realistic_in_docker.sh`
(R `n_core=4` vs Rust `--threads 8`) before Terra submission. See [UPSTREAM.md](UPSTREAM.md).

## CLI

```bash
tractor-mix-score \
  --null-export null_export \
  --dosage-files anc_00.dosage.txt.gz ... \
  --out phenotype.tractor_mix.tsv \
  --ac-threshold 50 \
  --threads 8 \
  --chunk-size 2048
```

Progress logs go to stderr (`variants/sec`, ETA). `null_export/` is written by
[`scripts/fit_null.R`](../../scripts/fit_null.R) after `glmmkin`.

## Output columns

`CHR POS ID REF ALT Chi2 P Eff_anc* SE_anc* Pval_anc* AC_count* include_anc*`

## Tests and benchmarks

```bash
cd tractor_mix/tractor_mix_score
cargo test
```

### Parity (exact match gate)

```bash
./scripts/run_oracle_parity_in_docker.sh   # requires tractor-mix-pilot:0.4.2 image
```

48-variant fixture: byte-identical TSV vs vendored R oracle.

### Realistic benchmark (pre-Terra gate)

```bash
./scripts/bench_realistic_in_docker.sh
```

Default: 4000 samples × 30k variants. **Must beat parallel R (`n_core=4`)** before production runs.

Env: `N_SAMPLES`, `N_SITES`, `R_CORES=4`, `RUST_THREADS=8`.

### Quick timing bench (toy scale)

```bash
./scripts/bench_r_vs_rust_in_docker.sh
```

## Docker

Built as a static musl binary in `tractor-mix-pilot:0.4.2` at `/usr/local/bin/tractor-mix-score`.

## Upstream contribution

See [UPSTREAM.md](UPSTREAM.md).
