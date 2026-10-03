# FELIX dosage-VCF and carrier-QC validation

Reproduces `felix/SV_SCORER_DESIGN.md` §3.1–3.2 on synthetic data (2,000
samples, two ancestries, binary 10% / 2% and quantitative traits). Rerun it
whenever `felix/patches/` is carried to a newer FELIX.

Runs under `linux/amd64` emulation on a Mac; expect roughly an hour end to end.
Use a work directory outside the repo (or under `scratch/`, which git ignores);
the VCFs are ~70 MB.

## Steps

```bash
REPO=$(git rev-parse --show-toplevel)
WORK=$REPO/scratch/dosage_qc && mkdir -p $WORK && cd $WORK
cp $REPO/felix/eval/dosage_qc/*.py $REPO/felix/eval/dosage_qc/*.sh .

# 1. Synthetic cohort, SNVs, repeat loci (needs numpy + pandas)
python3 simulate.py && python3 sim_null_rep.py

# 2. Stock FELIX baseline (nulls, FELIXla vs dosage VCF, 2,000 null repeat loci) -> out/
STOCK=us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.1.0
docker run --rm --platform linux/amd64 -v "$PWD":/w --entrypoint bash $STOCK /w/run.sh
docker run --rm --platform linux/amd64 -v "$PWD":/w --entrypoint bash $STOCK /w/run_null.sh
python3 compare.py      # §3.1: SNV equivalence, scale invariance, dropped tests
python3 calib.py        # §3.1: null calibration of the stock dosage path

# 3. Patched FELIX (felix/patches) -> val/
docker build --platform linux/amd64 -f $REPO/felix/eval/dosage_qc/Dockerfile.patch-test \
  -t felix-dosageqc-test:local $REPO
docker run --rm --platform linux/amd64 -v "$PWD":/w --entrypoint bash felix-dosageqc-test:local /w/validate.sh
python3 val_compare.py  # §3.2: switch off identical; switch on keeps all loci, calibrated

# 4. run_felix_step2.R command construction (stub step 2)
docker run --rm --platform linux/amd64 -v "$PWD":/w -v "$REPO":/repo --entrypoint bash \
  felix-dosageqc-test:local /w/validate_wrapper_stub.sh
```

`run.sh` / `run_null.sh` need `bcftools` and `plink2`, which the felix-pilot image
has and the bare FELIX image does not. `felix-pilot:0.1.0` is the last image
without the patch.

## Pass criteria

| Check | Pass |
|---|---|
| `val_compare.py` V1 | All three files `IDENTICAL` (switch off = stock FELIX) |
| `val_compare.py` V2 | Same SNVs tested; max \|Δlog10 p\| = 0 |
| `val_compare.py` V3 | No dropped tests at `minMAC` = 20 carriers; λ near 1; type I error near nominal |
| `validate_wrapper_stub.sh` | `run1 ok`, `run2 ok`, both refusal messages, `ENV carrier` only in VCF mode |
