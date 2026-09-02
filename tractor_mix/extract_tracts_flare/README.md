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

Annotate repeat-mediated SVs first, then extract (same CLI; `RU_TEST` sites
switch to copy-number dosage automatically):

```bash
python3 scripts/annotate_repeat_units.py \
  --vcf JOINT.vcf.gz \
  --simple-repeat-bed simpleRepeat.bed.gz \
  --out JOINT.ru.vcf.gz

mkdir -p OUTDIR
extract-tracts-flare \
  --vcf JOINT.ru.vcf.gz \
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
- `{prefix}.collapse.anc{i}.dosage.txt[.gz]` — biallelic 0/1 collapse of `RU_TEST` sites
- `{prefix}.split.anc{i}.dosage.txt[.gz]` — one 0/1 row per ALT of each `RU_TEST` site (`ID:1`, `ID:2`, …)
- `dosage_sample_order.txt` in the output directory

Sites with INFO `RU_TEST` use repeat-unit copy-number dosage: `C = CN_REF` on the
REF haplotype, or `C = CN_REF + sign(SVTYPE)×RU[a]` on ALT `a` (`sign` is −1 for
DEL and +1 for INS/DUP). Missing `CN_REF` is extra-units-vs-REF (`CN_REF=0`).
`dosage[anc] += C`. Multi-allelic VNTR-like records stay one locus in the main
dosage file. Non-`RU_TEST` sites keep the classic biallelic allele-`1` 0/1/2
dosage. Comparison encodings are written from the same haplotypes in one pass.

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
