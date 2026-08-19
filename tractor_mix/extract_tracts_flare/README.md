# extract-tracts-flare

Streaming Rust replacement for Atkinson-Lab Tractor
`extract_tracts_flare.py` (pinned SHA `2860728d`). It reads a FLARE
`GT:AN1:AN2` VCF and writes per-ancestry dosage and hapcount TSVs with the
same semantics as the Python script.

## CLI

```bash
extract-tracts-flare \
  --vcf INPUT.vcf.gz \
  --num-ancs 5 \
  --output-dir OUTDIR \
  --compress-output
```

Flags matching the Python script:

| Flag | Meaning |
|------|---------|
| `--vcf` | `.vcf` or `.vcf.gz` |
| `--num-ancs` | Ancestry count (`0 .. N-1`) |
| `--output-dir` | Must already exist (default: directory of the VCF) |
| `--output-vcf` | Also write ancestry-specific VCFs |
| `--compress-output` | zlib gzip (not bgzip) |

Additive flags:

| Flag | Meaning |
|------|---------|
| `--samples` / `--keep` | Subset and reorder sample columns (one ID per line) |
| `--threads` | Reserved; extract is currently single-threaded |

Outputs (prefix = VCF basename without `.vcf` / `.vcf.gz`):

- `{prefix}.anc{i}.dosage.txt[.gz]`
- `{prefix}.anc{i}.hapcount.txt[.gz]`
- `dosage_sample_order.txt` in the output directory

## Tests

```bash
cd tractor_mix/extract_tracts_flare
cargo test
```

Parity tests run the vendored Python oracle in `oracle/extract_tracts_flare.py`
against synthetic FLARE-like VCFs in `testdata/` and compare **decompressed**
text (gzip headers/mtime are ignored).

Generate a larger synthetic VCF:

```bash
python3 testdata/gen_flare_vcf.py --out /tmp/bench.vcf.gz \
  --sites 1000 --samples 1000 --num-ancs 5
cargo test -- --ignored   # 1k×1k smoke
```

Time Rust vs Python on the same file:

```bash
python3 testdata/gen_flare_vcf.py --out /tmp/bench.vcf.gz --sites 2000 --samples 2000
mkdir -p /tmp/py /tmp/rs
time python3 oracle/extract_tracts_flare.py --vcf /tmp/bench.vcf.gz --num-ancs 5 --output-dir /tmp/py --compress-output
time cargo run --release -- --vcf /tmp/bench.vcf.gz --num-ancs 5 --output-dir /tmp/rs --compress-output
```

Optional real-FORMAT smoke (not in CI): run FLARE’s
[test/run.flare.test](https://github.com/browning-lab/flare/blob/master/test/run.flare.test)
then point both extractors at `flare.out.anc.vcf.gz`. Do not commit AoU VCFs;
a chr22 slice is a Workbench-only capacity check.

## Docker

The Tractor-Mix image builds a **static musl** binary on Alpine (so it runs on
the older glibc in `wzhou88/saige:1.3.3`) and installs
`/usr/local/bin/extract-tracts-flare`. Image build runs `--help` on the SAIGE
base so a non-executable binary fails the build. `TractorMixPilot.wdl` calls it
instead of `/opt/Tractor/scripts/extract_tracts_flare.py`.
