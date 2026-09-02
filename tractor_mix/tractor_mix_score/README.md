# tractor-mix-score

Streaming Rust scorer for Tractor-Mix (GMMAT binomial, sparse GRM) and FELIX/SAIGE
Step 2 **normal-approximation** tests (quantitative + binomial; homogeneous /
heterogeneous / Cauchy combination). Dosages are parsed as `f64` (integer 0/1/2
files still work).

**Legacy GMMAT path:** same results as the pinned R oracle on the parity fixture.
`--threads` parallelizes variant scoring within each chunk (deterministic output
order). See [UPSTREAM.md](UPSTREAM.md).

Binary SPA is out of scope; FELIX mode uses the score-test chi-square (normal
approx) plus one scalar variance ratio.

## CLI

```bash
# GMMAT Tractor-Mix (auto-detected from meta.json source=gmmat / family=binomial)
tractor-mix-score \
  --null-export null_export \
  --dosage-files anc_00.dosage.txt.gz ... \
  --out phenotype.tractor_mix.tsv \
  --ac-threshold 50 \
  --mode auto \
  --threads 8 \
  --chunk-size 2048

# FELIX / SAIGE Step 1 export (f64 dosages, hom/het/CCT, no AC>50)
tractor-mix-score \
  --null-export null_export \
  --dosage-files anc_00.dosage.txt.gz ... \
  --out phenotype.felix.tsv \
  --mode felix \
  --variance-ratio 1.05 \
  --min-copy-var 1e-6 \
  --threads 8
```

`--mode auto` (default) uses FELIX scoring when `meta.json` has `source=felix`
(or `saige`) or `family` is `gaussian`/`quantitative`; otherwise the GMMAT
legacy path. `--ac-threshold` applies only in legacy mode. FELIX mode drops
ancestries with no dosage variance and applies one scalar variance ratio
(from `meta.json` / `variance_ratio.txt`, overridable with `--variance-ratio`).

Export a FELIX/SAIGE Step 1 `.rda` with
[`scripts/export_felix_null.R`](../../scripts/export_felix_null.R) or
`fit_null.R --step1-rda ...`. GMMAT exports still come from
[`scripts/fit_null.R`](../../scripts/fit_null.R) after `glmmkin`.

## Output columns

**Legacy (GMMAT):** `CHR POS ID REF ALT Chi2 P Eff_anc* SE_anc* Pval_anc* AC_count* include_anc*`

**FELIX:** `CHR POS MarkerID Allele1 Allele2` plus per-ancestry and `ancALL`
`AC_Allele2_* AF_Allele2_* BETA_* SE_* Tstat_* var_* p.value_*`, then
`P_het_admixed P_hom_admixed P_cct_admixed` and identical `_c` copies
(conditional analysis is not implemented; `_c` matches the unconditional tests).

## Tests and benchmarks

```bash
cd tractor_mix/tractor_mix_score
cargo test
```

### Parity (exact match gate, GMMAT)

```bash
./scripts/run_oracle_parity_in_docker.sh   # requires tractor-mix-pilot:0.4.2 image
```

48-variant fixture: byte-identical TSV vs vendored R oracle (`--mode legacy`).

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
