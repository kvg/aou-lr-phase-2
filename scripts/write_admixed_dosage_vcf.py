#!/usr/bin/env python3
"""Write FELIX admixed dosage VCF (DS{k}, ANC{k}) for repeat-mediated SVs.

Input: the joint phased VCF with FORMAT GT:AN1:AN2 (FLARE ancestry propagated
onto every site) and INFO from ``annotate_repeat_units.py`` (``RU_TEST``,
``RU`` per ALT, optional ``SVTYPE``). Only ``RU_TEST`` records are written.

Per haplotype ``h`` with allele ``a`` and ancestry ``k``::

    x_h = 0                     if a == 0 (REF)
    x_h = sign(SVTYPE) * RU[a]  otherwise (repeat units relative to REF)

    DS{k+1} = sum_{h: anc_h = k} x_h / s      ANC{k+1} = number of ancestry-k haplotypes

``s = max |x_h|`` at the locus (written as INFO ``RU_SCALE``), so
``|DS{k}| <= ANC{k} <= 2``. Score tests are invariant to ``s``; multiply FELIX
``BETA`` by ``1/s`` for an effect per repeat unit. The coding is never shifted:
a shift adds ``c * ANC{k}`` and turns the per-ancestry tests into local-ancestry
tests (felix/SV_SCORER_DESIGN.md).

Matches ``extract-tracts-flare --ru-baseline ref`` (same sign rule; ``CN_REF``
cancels). Run FELIX step 2 on the output with ``--vcfField=DS`` and
``FELIX_DOSAGE_QC=carrier`` (felix/patches/0001-dosage-carrier-qc.patch).

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


def allele_units(info: dict[str, Optional[str]], alt: str) -> list[float]:
    """x for allele index 0..n_alt (REF = 0)."""
    raw = info.get("RU")
    if not raw:
        raise ValueError("RU_TEST record without RU")
    ru = [float(v) for v in raw.split(",")]
    n_alt = alt.count(",") + 1
    if len(ru) != n_alt:
        raise ValueError(f"RU has {len(ru)} values for {n_alt} ALT alleles")
    sign = svtype_sign(info, alt)
    return [0.0] + [sign * v for v in ru]


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
) -> tuple[list[list[float]], list[list[int]], float]:
    """Per-sample (x sums per ancestry, haplotype counts per ancestry, max |x|)."""
    try:
        gi, a1i, a2i = fmt_keys.index("GT"), fmt_keys.index("AN1"), fmt_keys.index("AN2")
    except ValueError as exc:
        raise ValueError(f"{where}: FORMAT needs GT, AN1, AN2 (have {':'.join(fmt_keys)})") from exc
    ds_rows: list[list[float]] = []
    anc_rows: list[list[int]] = []
    max_abs = 0.0
    for field in sample_fields:
        parts = field.split(":")
        gt = parts[gi]
        if "|" not in gt:
            raise ValueError(f"{where}: unphased or missing GT {gt!r}")
        alleles = gt.split("|")
        ancs = (parts[a1i], parts[a2i])
        ds = [0.0] * num_ancs
        cnt = [0] * num_ancs
        for allele_s, anc_s in zip(alleles, ancs):
            if allele_s == "." or anc_s in (".", ""):
                raise ValueError(f"{where}: missing allele or ancestry in {field!r}")
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
    return ds_rows, anc_rows, max_abs


def header_lines(contigs: Sequence[str], num_ancs: int, samples: Sequence[str]) -> list[str]:
    out = ["##fileformat=VCFv4.2"]
    out.extend(contigs)
    out.append(
        '##INFO=<ID=RU_SCALE,Number=1,Type=Float,Description="DS{k} = sum of repeat units '
        'relative to REF over ancestry-k haplotypes / RU_SCALE">'
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
) -> dict[str, int]:
    stats = {"records_in": 0, "ru_records": 0, "written": 0, "monomorphic_skipped": 0}
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
                    missing = [s for s in keep if s not in idx]
                    if missing:
                        raise SystemExit(f"{len(missing)} keep-list samples not in VCF, e.g. {missing[:3]}")
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
            units = allele_units(info, f[4])
            fields = [f[9 + i] for i in order]
            ds_rows, anc_rows, max_abs = site_dosages(fields, f[8].split(":"), units, num_ancs, where)
            if max_abs == 0.0:
                stats["monomorphic_skipped"] += 1
                continue
            cells = []
            for ds, cnt in zip(ds_rows, anc_rows):
                cells.append(":".join([fmt_num(v / max_abs) for v in ds] + [str(c) for c in cnt]))
            info_out = f"RU_SCALE={fmt_num(max_abs)}"
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
    args = p.parse_args(argv)

    keep = None
    if args.samples:
        keep = [ln.strip() for ln in args.samples.read_text().splitlines() if ln.strip()]
    out = sys.stdout if args.out == "-" else open(args.out, "w")
    try:
        stats = write_admixed_vcf(args.vcf, out, num_ancs=args.num_ancs, keep=keep)
    finally:
        if out is not sys.stdout:
            out.close()
    print(json.dumps(stats), file=sys.stderr)
    if args.stats_json:
        args.stats_json.write_text(json.dumps(stats, indent=2) + "\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
