#!/usr/bin/env python3
"""Build an include_sites BED of HiFi↔srWGS genotype-concordant variants.

Compares two VCFs on the intersection of samples (and optional keep-list),
restricted to biallelic SNVs in ``--region``. A site is kept when the fraction
of sample pairs with non-missing genotypes on both sides that match is at least
``--min-concordance`` (default 1.0 = perfect agreement among callable pairs),
and at least ``--min-called`` samples are callable on both sides.

Writes:
  * ``*.bed.gz`` (+ ``.tbi`` when bgzip/tabix are available) for FLARE
    ``include_sites`` / bcftools ``-T``
  * ``*.tsv.gz`` with chrom/pos/ref/alt + concordance stats (provenance)

Example::

    python3 scripts/flare_build_concordance_sites.py \\
      --hifi-vcf aou_lr_phase2_v1.chr20.vcf.gz \\
      --sr-vcf srwgs.chr20.vcf.gz \\
      --region chr20 \\
      --out-prefix hifi_srwgs_concordant.chr20
"""

from __future__ import annotations

import argparse
import gzip
import shutil
import subprocess
import sys
from pathlib import Path


def run(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    print("+", " ".join(map(str, cmd)), file=sys.stderr)
    return subprocess.run(cmd, check=True, text=True, **kwargs)


def bcftools_samples(vcf: str) -> list[str]:
    out = run(["bcftools", "query", "-l", vcf], capture_output=True)
    return [s for s in out.stdout.splitlines() if s.strip()]


def load_keep(path: Path | None) -> set[str] | None:
    if path is None:
        return None
    ids: set[str] = set()
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            ids.add(line.split()[0])
    return ids


def open_gt_stream(
    vcf: str,
    *,
    region: str,
    samples: list[str],
    sample_file: Path,
) -> subprocess.Popen[str]:
    sample_file.write_text("\n".join(samples) + "\n")
    cmd = [
        "bcftools",
        "query",
        "-S",
        str(sample_file),
        "-m2",
        "-M2",
        "-v",
        "snps",
        "-f",
        "%CHROM\t%POS\t%REF\t%ALT[\t%GT]\n",
    ]
    if region.strip():
        cmd.extend(["-r", region.strip()])
    cmd.append(vcf)
    print("+", " ".join(cmd), file=sys.stderr)
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def parse_row(line: str) -> tuple[str, int, str, str, list[str]]:
    parts = line.rstrip("\n").split("\t")
    if len(parts) < 5:
        raise ValueError(f"short query row: {line[:80]!r}")
    chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
    gts = parts[4:]
    return chrom, int(pos_s), ref, alt, gts


def norm_gt(gt: str) -> str | None:
    if gt in {".", "./.", ".|.", ""}:
        return None
    alleles = gt.replace("|", "/").split("/")
    if len(alleles) != 2 or any(a == "." for a in alleles):
        return None
    try:
        a, b = int(alleles[0]), int(alleles[1])
    except ValueError:
        return None
    if a > b:
        a, b = b, a
    return f"{a}/{b}"


def chrom_key(chrom: str) -> tuple:
    c = chrom[3:] if chrom.startswith("chr") else chrom
    if c.isdigit():
        return (0, int(c))
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}
    return (1, order.get(c.upper(), 100), c)


def site_key(chrom: str, pos: int, ref: str, alt: str) -> tuple:
    return (chrom_key(chrom), pos, ref, alt)


def compress_bed(sorted_bed: Path, out_gz: Path) -> None:
    if out_gz.exists():
        out_gz.unlink()
    tbi = Path(str(out_gz) + ".tbi")
    if tbi.exists():
        tbi.unlink()
    if shutil.which("bgzip"):
        with out_gz.open("wb") as fh:
            subprocess.run(["bgzip", "-c", str(sorted_bed)], stdout=fh, check=True)
        if shutil.which("tabix"):
            run(["tabix", "-p", "bed", str(out_gz)])
    else:
        with sorted_bed.open("rb") as src, gzip.open(out_gz, "wb") as dest:
            shutil.copyfileobj(src, dest)


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hifi-vcf", required=True, help="Long-read / HiFi target VCF")
    p.add_argument("--sr-vcf", required=True, help="Short-read WGS VCF for the same cohort")
    p.add_argument("--region", default="", help="chr or chr:start-end (empty = whole files)")
    p.add_argument("--samples", type=Path, default=None, help="Optional sample keep-list")
    p.add_argument("--min-concordance", type=float, default=1.0)
    p.add_argument("--min-called", type=int, default=10)
    p.add_argument("--out-prefix", type=Path, required=True)
    args = p.parse_args()

    if not (0.0 < args.min_concordance <= 1.0):
        raise SystemExit("--min-concordance must be in (0, 1]")
    if args.min_called < 1:
        raise SystemExit("--min-called must be >= 1")

    hifi_samples = set(bcftools_samples(args.hifi_vcf))
    sr_samples = set(bcftools_samples(args.sr_vcf))
    shared = sorted(hifi_samples & sr_samples)
    keep = load_keep(args.samples)
    if keep is not None:
        shared = [s for s in shared if s in keep]
    if len(shared) < args.min_called:
        raise SystemExit(
            f"only {len(shared)} shared samples after filters; "
            f"need >= {args.min_called}"
        )
    print(f"shared samples: {len(shared)}", file=sys.stderr)

    args.out_prefix.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = args.out_prefix.parent / f".{args.out_prefix.name}.tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    hifi_sf = tmp_dir / "hifi.samples"
    sr_sf = tmp_dir / "sr.samples"

    hifi_proc = open_gt_stream(
        args.hifi_vcf, region=args.region, samples=shared, sample_file=hifi_sf
    )
    sr_proc = open_gt_stream(
        args.sr_vcf, region=args.region, samples=shared, sample_file=sr_sf
    )
    assert hifi_proc.stdout is not None and sr_proc.stdout is not None

    bed_path = tmp_dir / "sites.bed"
    tsv_path = Path(str(args.out_prefix) + ".tsv.gz")
    n_hifi = n_sr = n_join = n_keep = 0
    hifi_line = hifi_proc.stdout.readline()
    sr_line = sr_proc.stdout.readline()
    hifi_row = parse_row(hifi_line) if hifi_line else None
    sr_row = parse_row(sr_line) if sr_line else None

    with bed_path.open("w") as bed_fh, gzip.open(tsv_path, "wt") as tsv_fh:
        tsv_fh.write(
            "chrom\tpos\tref\talt\tn_called\tn_concordant\tconcordance\tn_samples\n"
        )
        while hifi_row is not None and sr_row is not None:
            hk = site_key(*hifi_row[:4])
            sk = site_key(*sr_row[:4])
            if hk < sk:
                n_hifi += 1
                hifi_line = hifi_proc.stdout.readline()
                hifi_row = parse_row(hifi_line) if hifi_line else None
                continue
            if sk < hk:
                n_sr += 1
                sr_line = sr_proc.stdout.readline()
                sr_row = parse_row(sr_line) if sr_line else None
                continue

            n_join += 1
            chrom, pos, ref, alt, gts_h = hifi_row
            _, _, _, _, gts_s = sr_row
            if len(gts_h) != len(shared) or len(gts_s) != len(shared):
                raise SystemExit("GT column count mismatch vs shared sample list")
            called = conc = 0
            for gh, gs in zip(gts_h, gts_s):
                nh, ns = norm_gt(gh), norm_gt(gs)
                if nh is None or ns is None:
                    continue
                called += 1
                if nh == ns:
                    conc += 1
            frac = (conc / called) if called else 0.0
            if called >= args.min_called and frac >= args.min_concordance:
                n_keep += 1
                bed_fh.write(f"{chrom}\t{pos - 1}\t{pos}\n")
                tsv_fh.write(
                    f"{chrom}\t{pos}\t{ref}\t{alt}\t{called}\t{conc}\t{frac:.6f}\t{len(shared)}\n"
                )

            n_hifi += 1
            n_sr += 1
            hifi_line = hifi_proc.stdout.readline()
            sr_line = sr_proc.stdout.readline()
            hifi_row = parse_row(hifi_line) if hifi_line else None
            sr_row = parse_row(sr_line) if sr_line else None

        while hifi_row is not None:
            n_hifi += 1
            hifi_line = hifi_proc.stdout.readline()
            hifi_row = parse_row(hifi_line) if hifi_line else None
        while sr_row is not None:
            n_sr += 1
            sr_line = sr_proc.stdout.readline()
            sr_row = parse_row(sr_line) if sr_line else None

    hifi_err = hifi_proc.communicate()[1]
    sr_err = sr_proc.communicate()[1]
    if hifi_proc.returncode:
        sys.stderr.write(hifi_err or "")
        raise SystemExit(f"bcftools query (hifi) failed: {hifi_proc.returncode}")
    if sr_proc.returncode:
        sys.stderr.write(sr_err or "")
        raise SystemExit(f"bcftools query (sr) failed: {sr_proc.returncode}")

    bed_gz = Path(str(args.out_prefix) + ".bed.gz")
    compress_bed(bed_path, bed_gz)

    print(
        f"hifi_sites_seen={n_hifi} sr_sites_seen={n_sr} "
        f"allele_matched={n_join} kept={n_keep}",
        file=sys.stderr,
    )
    print(f"Wrote {bed_gz}", file=sys.stderr)
    print(f"Wrote {tsv_path}", file=sys.stderr)
    if n_keep == 0:
        raise SystemExit("no concordant sites retained; check inputs / thresholds")

    # Best-effort temp cleanup
    for path in (bed_path, hifi_sf, sr_sf):
        path.unlink(missing_ok=True)
    try:
        tmp_dir.rmdir()
    except OSError:
        pass


if __name__ == "__main__":
    main()
