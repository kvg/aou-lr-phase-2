#!/usr/bin/env python3
"""Write FELIX admixed dosage VCF (DS{k}, ANC{k}) for repeat-mediated SVs.

Input: the joint phased VCF with FORMAT GT:AN1:AN2 (FLARE ancestry propagated
onto every site) and INFO from ``annotate_repeat_units.py`` (``RU_TEST``,
``RU`` per ALT, optional ``SVTYPE``) or from ``aggregate_repeat_loci.py``
(``RU_TEST``, signed ``RU_DOSAGE`` per ALT). Only ``RU_TEST`` records are written.

Per haplotype ``h`` with allele ``a`` and ancestry ``k``::

    x_h = 0                     if a == 0 (REF)
    x_h = RU_DOSAGE[a]          when the record has RU_DOSAGE (one locus per record)
    x_h = sign(a) * RU[a]       otherwise (repeat units relative to REF)

``sign(a)`` is -1 when the ALT is shorter than REF and +1 when longer; symbolic
ALTs fall back to ``SVTYPE`` (DEL = -1).

    DS{k+1} = sum_{h: anc_h = k} x_h / s      ANC{k+1} = number of ancestry-k haplotypes

``s = max |x_h|`` at the locus (written as INFO ``RU_SCALE``), so
``|DS{k}| <= ANC{k} <= 2``. Score tests are invariant to ``s``; multiply FELIX
``BETA`` by ``1/s`` for an effect per repeat unit. The coding is never shifted:
a shift adds ``c * ANC{k}`` and turns the per-ancestry tests into local-ancestry
tests (felix/SV_SCORER_DESIGN.md).

Matches ``extract-tracts-flare --ru-baseline ref`` on ``annotate_repeat_units.py``
output, which writes ``SVTYPE`` for every ``RU_TEST`` record (``CN_REF`` cancels). Run FELIX step 2 on the output with ``--vcfField=DS`` and
``FELIX_DOSAGE_QC=carrier`` (felix/patches/0001-dosage-carrier-qc.patch).

Missing haplotypes
------------------
``aggregate_repeat_loci.py`` writes ``.|.`` for a haplotype whose contributing
record was missing and ``./.`` for a sample whose length changes cannot be
assigned to haplotypes, so real loci arrive with missing calls.
``--missing ref`` (default) gives such a haplotype ``x = 0`` and **excludes it
from** ``ANC{k}``; ``--max-hap-missing`` drops a locus above a missingness
fraction, recorded per locus as INFO ``RU_HAP_MISS``. ``--missing error`` keeps
the original strict refusal.

Excluding the haplotype from ``ANC{k}`` rather than counting it is what keeps
the FELIX carrier-QC patch exact. That patch counts a sample as a carrier when
``DS{k} != 0``, among samples with ``ANC{k} > 0``; a missing haplotype is
neither a carrier nor part of the denominator, which is the same convention
``aggregate_repeat_loci.py`` uses for its ``n_hap`` tallies. Note that the patch
evaluates carrier counts at two points in ``mainMarkerAdmixedInCPP`` — once
before FELIX's own imputation and once after
(``TRACTORHYBRID_TIMED_IMPUTE_FAKEFLIP``). Emitting ``0`` rather than a missing
``DS`` means both evaluations see the same value and FELIX's ``g_impute_method``
never fires on these records, so the two call sites cannot disagree.

Example::

    python3 scripts/write_admixed_dosage_vcf.py \\
      --vcf chr22.gt_an.ru.vcf.gz --num-ancs 5 --samples analysis_samples.txt \\
      --out chr22.ru.admixed.vcf
    bgzip chr22.ru.admixed.vcf && bcftools index -c chr22.ru.admixed.vcf.gz
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from pathlib import Path
from typing import Optional, Sequence, TextIO


def open_text(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def info_dict(info: str) -> dict[str, Optional[str]]:
    out: dict[str, Optional[str]] = {}
    if info in ("", "."):
        return out
    for item in info.split(";"):
        key, sep, val = item.partition("=")
        out[key] = val if sep else None
    return out


def is_ru_test(info: dict[str, Optional[str]]) -> bool:
    if "RU_TEST" not in info:
        return False
    return info["RU_TEST"] not in ("0", "False", "false")


def svtype_sign(info: dict[str, Optional[str]], alt: str) -> int:
    """-1 for deletions, +1 otherwise (same rule as extract-tracts-flare)."""
    svtype = info.get("SVTYPE")
    if svtype is not None:
        first = svtype.replace(":", ",").split(",")[0]
        return -1 if first.upper() == "DEL" else 1
    return -1 if alt[:4].upper() == "<DEL" else 1


def allele_sign(info: dict[str, Optional[str]], ref: str, alt: str) -> int:
    if not alt.startswith("<") and "[" not in alt and "]" not in alt and alt != "*":
        if len(alt) < len(ref):
            return -1
        if len(alt) > len(ref):
            return 1
    return svtype_sign(info, alt)


def allele_units(info: dict[str, Optional[str]], ref: str, alt: str) -> list[float]:
    """x for allele index 0..n_alt (REF = 0)."""
    alts = alt.split(",")
    signed = info.get("RU_DOSAGE")
    if signed:
        vals = [float(v) for v in signed.split(",")]
        if len(vals) != len(alts):
            raise ValueError(f"RU_DOSAGE has {len(vals)} values for {len(alts)} ALT alleles")
        return [0.0] + vals
    raw = info.get("RU")
    if not raw:
        raise ValueError("RU_TEST record without RU or RU_DOSAGE")
    ru = [float(v) for v in raw.split(",")]
    if len(ru) != len(alts):
        raise ValueError(f"RU has {len(ru)} values for {len(alts)} ALT alleles")
    return [0.0] + [allele_sign(info, ref, a) * v for a, v in zip(alts, ru)]


def fmt_num(v: float) -> str:
    if v == int(v):
        return str(int(v))
    return f"{v:.8g}"


def site_dosages(
    sample_fields: Sequence[str],
    fmt_keys: Sequence[str],
    units: Sequence[float],
    num_ancs: int,
    where: str,
    missing: str = "ref",
) -> tuple[list[list[float]], list[list[int]], float, int, int]:
    """Per-sample (x sums per ancestry, haplotype counts per ancestry, max |x|,
    missing haplotypes, total haplotypes).

    A haplotype is missing when its allele or its ancestry is ``.``, or when the
    genotype is unphased (``/``) — the locus aggregator writes ``./.`` for a
    sample whose length changes cannot be assigned to haplotypes, and an
    unphased heterozygote has no haplotype-to-ancestry mapping to test.

    ``missing``:
      ``error`` — refuse the record (the original strict behaviour).
      ``ref``   — the haplotype contributes ``x = 0`` and is **excluded** from
                  ``ANC{k}``. Under REF-relative coding that keeps carrier QC
                  exact: a carrier is a sample with a non-zero ``DS{k}``, i.e.
                  one that carries a non-reference repeat allele, and a missing
                  haplotype is neither a carrier nor part of the denominator.
    """
    if missing not in ("error", "ref"):
        raise ValueError(f"missing must be 'error' or 'ref', got {missing!r}")
    try:
        gi, a1i, a2i = fmt_keys.index("GT"), fmt_keys.index("AN1"), fmt_keys.index("AN2")
    except ValueError as exc:
        raise ValueError(f"{where}: FORMAT needs GT, AN1, AN2 (have {':'.join(fmt_keys)})") from exc
    ds_rows: list[list[float]] = []
    anc_rows: list[list[int]] = []
    max_abs = 0.0
    n_missing = 0
    n_hap = 0
    for field in sample_fields:
        parts = field.split(":")
        gt = parts[gi]
        phased = "|" in gt
        if not phased and missing == "error":
            raise ValueError(f"{where}: unphased or missing GT {gt!r}")
        alleles = gt.replace("/", "|").split("|")
        ancs = (parts[a1i], parts[a2i])
        ds = [0.0] * num_ancs
        cnt = [0] * num_ancs
        for allele_s, anc_s in zip(alleles, ancs):
            n_hap += 1
            if allele_s in (".", "") or anc_s in (".", "") or not phased:
                if missing == "error":
                    raise ValueError(f"{where}: missing allele or ancestry in {field!r}")
                n_missing += 1
                continue
            allele, anc = int(allele_s), int(anc_s)
            if not 0 <= anc < num_ancs:
                raise ValueError(f"{where}: ancestry {anc} outside 0..{num_ancs - 1}")
            if allele >= len(units):
                raise ValueError(f"{where}: allele {allele} has no RU value")
            x = units[allele]
            ds[anc] += x
            cnt[anc] += 1
            max_abs = max(max_abs, abs(x))
        ds_rows.append(ds)
        anc_rows.append(cnt)
    return ds_rows, anc_rows, max_abs, n_missing, n_hap


def header_lines(contigs: Sequence[str], num_ancs: int, samples: Sequence[str]) -> list[str]:
    out = ["##fileformat=VCFv4.2"]
    out.extend(contigs)
    out.append(
        '##INFO=<ID=RU_SCALE,Number=1,Type=Float,Description="DS{k} = sum of repeat units '
        'relative to REF over ancestry-k haplotypes / RU_SCALE">'
    )
    out.append(
        '##INFO=<ID=RU_HAP_MISS,Number=1,Type=Float,Description="Fraction of haplotypes with a '
        'missing allele or ancestry, or an unphased genotype; they contribute 0 to DS{k} and are '
        'excluded from ANC{k}">'
    )
    for k in range(1, num_ancs + 1):
        out.append(
            f'##FORMAT=<ID=DS{k},Number=1,Type=Float,Description="Scaled repeat-unit dosage, ancestry {k}">'
        )
    for k in range(1, num_ancs + 1):
        out.append(
            f'##FORMAT=<ID=ANC{k},Number=1,Type=Float,Description="Number of ancestry-{k} haplotypes">'
        )
    out.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    return out


def write_admixed_vcf(
    vcf: str,
    out: TextIO,
    *,
    num_ancs: int,
    keep: Optional[list[str]] = None,
    missing: str = "ref",
    max_hap_missing: float = 0.1,
) -> dict[str, int]:
    stats = {
        "records_in": 0, "ru_records": 0, "written": 0, "monomorphic_skipped": 0,
        "missingness_skipped": 0, "hap_missing": 0, "hap_total": 0,
    }
    contigs: list[str] = []
    order: Optional[list[int]] = None
    fmt_out = ":".join([f"DS{k}" for k in range(1, num_ancs + 1)] + [f"ANC{k}" for k in range(1, num_ancs + 1)])
    with open_text(vcf) as fh:
        for line in fh:
            if line.startswith("##"):
                if line.startswith("##contig="):
                    contigs.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                vcf_samples = cols[9:]
                if keep is None:
                    order = list(range(len(vcf_samples)))
                    names = vcf_samples
                else:
                    idx = {s: i for i, s in enumerate(vcf_samples)}
                    absent = [s for s in keep if s not in idx]
                    if absent:
                        raise SystemExit(f"{len(absent)} keep-list samples not in VCF, e.g. {absent[:3]}")
                    order = [idx[s] for s in keep]
                    names = keep
                out.write("\n".join(header_lines(contigs, num_ancs, names)) + "\n")
                continue
            if order is None:
                raise SystemExit("data before #CHROM header")
            stats["records_in"] += 1
            f = line.rstrip("\n").split("\t")
            info = info_dict(f[7])
            if not is_ru_test(info):
                continue
            stats["ru_records"] += 1
            where = f"{f[0]}:{f[1]}"
            units = allele_units(info, f[3], f[4])
            fields = [f[9 + i] for i in order]
            ds_rows, anc_rows, max_abs, n_miss, n_hap = site_dosages(
                fields, f[8].split(":"), units, num_ancs, where, missing=missing
            )
            if max_abs == 0.0:
                stats["monomorphic_skipped"] += 1
                continue
            miss_frac = (n_miss / n_hap) if n_hap else 0.0
            if miss_frac > max_hap_missing:
                stats["missingness_skipped"] += 1
                continue
            stats["hap_missing"] += n_miss
            stats["hap_total"] += n_hap
            cells = []
            for ds, cnt in zip(ds_rows, anc_rows):
                cells.append(":".join([fmt_num(v / max_abs) for v in ds] + [str(c) for c in cnt]))
            info_out = f"RU_SCALE={fmt_num(max_abs)};RU_HAP_MISS={fmt_num(round(miss_frac, 6))}"
            out.write("\t".join([f[0], f[1], f[2], f[3], f[4], ".", f[6], info_out, fmt_out, *cells]) + "\n")
            stats["written"] += 1
    return stats


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="Joint GT:AN1:AN2 VCF (.vcf/.vcf.gz, or - for stdin)")
    p.add_argument("--num-ancs", required=True, type=int)
    p.add_argument("--samples", type=Path, default=None, help="Keep-list; output columns in this order")
    p.add_argument("--out", required=True, help="Output VCF (plain text; bgzip afterwards), or - for stdout")
    p.add_argument("--stats-json", type=Path, default=None)
    p.add_argument(
        "--missing", choices=("ref", "error"), default="ref",
        help="ref: a missing or unphased haplotype contributes 0 to DS{k} and is excluded from "
             "ANC{k}; error: refuse any such record",
    )
    p.add_argument(
        "--max-hap-missing", type=float, default=0.1,
        help="Skip a locus whose missing-haplotype fraction exceeds this (default 0.1)",
    )
    args = p.parse_args(argv)
    if not 0.0 <= args.max_hap_missing <= 1.0:
        raise SystemExit("--max-hap-missing must be between 0 and 1")

    keep = None
    if args.samples:
        keep = [ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()]
    out = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        stats = write_admixed_vcf(
            args.vcf, out, num_ancs=args.num_ancs, keep=keep,
            missing=args.missing, max_hap_missing=args.max_hap_missing,
        )
    finally:
        if out is not sys.stdout:
            out.close()
    print(json.dumps(stats), file=sys.stderr)
    if args.stats_json:
        args.stats_json.write_text(json.dumps(stats, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
