#!/usr/bin/env python3
"""Build an include_sites BED of high-quality DeepVariant/GLnexus SNVs.

Keeps biallelic SNVs where the fraction of samples failing call QC is at most
``--max-frac-fail`` (default 0.05). A sample fails when any of:

* ``GQ < --min-gq`` (default 20), when GQ is present
* ``DP < --min-dp`` (default 10), when DP is present
* ``RNC`` contains ``I`` (incomplete gVCF / no-call), when ``--drop-rnc-i``

Sites missing both GQ and DP for every sample are dropped (cannot QC).

Writes:
  * ``*.bed.gz`` (+ ``.tbi`` when bgzip/tabix are available) for FLARE
    ``include_sites`` / bcftools ``-T``
  * ``*.tsv.gz`` with chrom/pos/ref/alt + fail fraction (provenance)

Example::

    python3 scripts/flare_build_call_qc_sites.py \\
      --vcf glnexus.chr22.vcf.gz \\
      --region chr22:26897597-36897597 \\
      --out-prefix call_qc.chr22_10mb
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


def parse_int(tok: str) -> int | None:
    tok = tok.strip()
    if not tok or tok in {".", "./.", ".|."}:
        return None
    try:
        return int(tok)
    except ValueError:
        return None


def rnc_has_i(rnc: str) -> bool:
    return "I" in (rnc or "").upper()


def sample_fails(
    gq_s: str,
    dp_s: str,
    rnc_s: str,
    *,
    min_gq: int,
    min_dp: int,
    drop_rnc_i: bool,
) -> bool | None:
    """Return True/False if QC-able, None if both GQ and DP missing."""
    gq = parse_int(gq_s)
    dp = parse_int(dp_s)
    if gq is None and dp is None and not (drop_rnc_i and rnc_has_i(rnc_s)):
        return None
    if gq is not None and gq < min_gq:
        return True
    if dp is not None and dp < min_dp:
        return True
    if drop_rnc_i and rnc_has_i(rnc_s):
        return True
    return False


def open_qc_stream(
    vcf: str,
    *,
    region: str,
    samples: list[str],
    sample_file: Path,
) -> subprocess.Popen[str]:
    sample_file.write_text("\n".join(samples) + "\n")
    # Per sample: GQ, DP, RNC
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
        "%CHROM\t%POS\t%REF\t%ALT[\t%GQ\t%DP\t%RNC]\n",
    ]
    if region.strip():
        cmd.extend(["-r", region.strip()])
    cmd.append(vcf)
    print("+", " ".join(cmd), file=sys.stderr)
    return subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        text=True,
        bufsize=1,
    )


def bgzip_tabix(bed: Path) -> Path:
    gz = Path(str(bed) + ".gz") if not str(bed).endswith(".gz") else bed
    if not str(bed).endswith(".gz"):
        run(["bgzip", "-f", str(bed)])
        bed = gz
    else:
        bed = gz
    if shutil.which("tabix"):
        run(["tabix", "-f", "-p", "bed", str(bed)])
    return bed


def main(argv: list[str] | None = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, help="DeepVariant/GLnexus joint VCF with GQ/DP/RNC")
    p.add_argument("--region", default="", help="chr or chr:start-end (optional)")
    p.add_argument("--keep", type=Path, help="Optional sample keep-list")
    p.add_argument("--min-gq", type=int, default=20)
    p.add_argument("--min-dp", type=int, default=10)
    p.add_argument(
        "--max-frac-fail",
        type=float,
        default=0.05,
        help="Keep site if fail_frac <= this (default 0.05)",
    )
    p.add_argument(
        "--drop-rnc-i",
        action=argparse.BooleanOptionalAction,
        default=True,
        help="Treat RNC containing I as a fail (default: true)",
    )
    p.add_argument("--out-prefix", required=True, help="Output path prefix")
    args = p.parse_args(argv)

    if not shutil.which("bcftools"):
        raise SystemExit("bcftools is required")

    samples = bcftools_samples(args.vcf)
    keep = load_keep(args.keep)
    if keep is not None:
        samples = [s for s in samples if s in keep]
    if not samples:
        raise SystemExit("no samples after keep-list filter")

    out_prefix = Path(args.out_prefix)
    out_prefix.parent.mkdir(parents=True, exist_ok=True)
    bed_path = Path(str(out_prefix) + ".bed")
    tsv_path = Path(str(out_prefix) + ".tsv.gz")
    sample_file = Path(str(out_prefix) + ".samples.txt")

    n_sites = 0
    n_kept = 0
    n_no_qc = 0
    proc = open_qc_stream(
        args.vcf, region=args.region, samples=samples, sample_file=sample_file
    )
    assert proc.stdout is not None

    with bed_path.open("wt") as bed_fh, gzip.open(tsv_path, "wt") as tsv_fh:
        tsv_fh.write(
            "chrom\tpos\tref\talt\tn_samples\tn_qcable\tn_fail\tfail_frac\tkeep\n"
        )
        for line in proc.stdout:
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 4:
                continue
            chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
            toks = parts[4:]
            if len(toks) != 3 * len(samples):
                continue
            n_sites += 1
            n_fail = 0
            n_qcable = 0
            for i in range(len(samples)):
                gq_s, dp_s, rnc_s = toks[3 * i], toks[3 * i + 1], toks[3 * i + 2]
                failed = sample_fails(
                    gq_s,
                    dp_s,
                    rnc_s,
                    min_gq=args.min_gq,
                    min_dp=args.min_dp,
                    drop_rnc_i=args.drop_rnc_i,
                )
                if failed is None:
                    continue
                n_qcable += 1
                if failed:
                    n_fail += 1
            if n_qcable == 0:
                n_no_qc += 1
                fail_frac = 1.0
                keep_site = False
            else:
                fail_frac = n_fail / n_qcable
                keep_site = fail_frac <= args.max_frac_fail
            tsv_fh.write(
                f"{chrom}\t{pos_s}\t{ref}\t{alt}\t{len(samples)}\t{n_qcable}\t"
                f"{n_fail}\t{fail_frac:.6f}\t{int(keep_site)}\n"
            )
            if keep_site:
                # BED half-open
                pos = int(pos_s)
                bed_fh.write(f"{chrom}\t{pos - 1}\t{pos}\n")
                n_kept += 1

    rc = proc.wait()
    if rc != 0:
        raise SystemExit(f"bcftools query failed with exit {rc}")

    bed_gz = bgzip_tabix(bed_path)
    print(
        f"sites={n_sites} kept={n_kept} no_qc={n_no_qc} "
        f"min_gq={args.min_gq} min_dp={args.min_dp} "
        f"max_frac_fail={args.max_frac_fail} -> {bed_gz}",
        file=sys.stderr,
    )
    print(f"stats -> {tsv_path}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
