#!/usr/bin/env python3
"""
Annotate a joint phased VCF with repeat-unit INFO tags.

SVs live in the same VCF as SNVs/indels (GT + AN1/AN2). No TRGT. Writes:

  RU_TEST  Flag      INS/DEL/DUP overlapping a tandem array with a usable period
  PERIOD   Integer   period (bp) from simpleRepeat/TRF or allele sequence
  MOTIF    String    consensus motif
  CN_REF   Float     reference copy number when known; omitted = extra-units-vs-REF
  RU       Integer A repeat units of each ALT (SVLEN/PERIOD for symbolic ALTs)

INV, BND, SNVs, and unique-sequence events are not RU_TEST.

``--simple-repeat-bed`` accepts:
  - 3-column BED (chrom/start/end), as staged by sv_01 / sv_annotation fixtures
  - UCSC simpleRepeat.txt (optional bin + period + copyNum + sequence)
  - TRF-style BED with period / copy-number / motif columns after chrom/start/end
"""

from __future__ import annotations

import argparse
import gzip
import math
import re
import sys
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Optional, TextIO

sys.path.insert(0, str(Path(__file__).resolve().parent))

from sv_site_utils import open_text, parse_info  # noqa: E402

RU_ELIGIBLE = frozenset({"INS", "DEL", "DUP"})
RU_EXCLUDED = frozenset({"INV", "BND", "SNV", "MNV", "CNV", "CPX", "CTX"})

INFO_HEADERS = [
    '##INFO=<ID=RU_TEST,Number=0,Type=Flag,Description="INS/DEL/DUP overlapping a tandem array with a usable period; eligible for repeat-unit dosage">',
    '##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Tandem repeat period in bp from simpleRepeat/TRF or allele sequence">',
    '##INFO=<ID=MOTIF,Number=1,Type=String,Description="Tandem repeat consensus motif">',
    '##INFO=<ID=CN_REF,Number=1,Type=Float,Description="Reference copy number of the overlapping tandem array; omitted means extra-units-vs-REF">',
    '##INFO=<ID=RU,Number=A,Type=Integer,Description="Repeat-unit count of each ALT (SVLEN/PERIOD for symbolic ALTs)">',
]

_CHROM_RE = re.compile(r"^(chr)?([1-9]|1[0-9]|2[0-2]|[XYM]|MT)$", re.IGNORECASE)
_DNA_RE = re.compile(r"^[ACGTNacgtn]+$")


@dataclass(frozen=True)
class TandemArray:
    start: int
    end: int
    period: Optional[int]
    motif: Optional[str]
    copy_num: Optional[float]


@dataclass(frozen=True)
class RuAnnotation:
    period: int
    motif: Optional[str]
    cn_ref: Optional[float]
    ru: tuple[int, ...]


def _is_chrom(token: str) -> bool:
    if not token:
        return False
    low = token.lower()
    if low.startswith("chr"):
        return True
    return bool(_CHROM_RE.match(token))


def _parse_int(token: str) -> Optional[int]:
    try:
        return int(token)
    except (TypeError, ValueError):
        return None


def _parse_float(token: str) -> Optional[float]:
    try:
        v = float(token)
    except (TypeError, ValueError):
        return None
    if math.isnan(v) or math.isinf(v):
        return None
    return v


def _looks_dna(token: str) -> bool:
    return bool(token) and bool(_DNA_RE.match(token))


def parse_tandem_bed_line(line: str) -> Optional[tuple[str, TandemArray]]:
    """Parse one simpleRepeat / TRF / 3-column BED line."""
    raw = line.rstrip("\n")
    if not raw or raw.startswith("#") or raw.startswith("track"):
        return None
    parts = raw.split("\t")
    if len(parts) == 1:
        parts = raw.split()
    if len(parts) < 3:
        return None

    offset = 0
    if not _is_chrom(parts[0]) and _parse_int(parts[0]) is not None and len(parts) >= 4:
        if _is_chrom(parts[1]):
            offset = 1
    chrom = parts[offset]
    start = _parse_int(parts[offset + 1])
    end = _parse_int(parts[offset + 2])
    if start is None or end is None:
        return None
    if end < start:
        start, end = end, start
    if end == start:
        end = start + 1

    rest = parts[offset + 3 :]
    period: Optional[int] = None
    copy_num: Optional[float] = None
    motif: Optional[str] = None

    if not rest:
        return chrom, TandemArray(start, end, None, None, None)

    # chrom start end period [copyNum] [motif/consensus...]
    if _parse_int(rest[0]) is not None and (
        len(rest) == 1 or _parse_float(rest[1]) is not None or _looks_dna(rest[1])
    ):
        period = _parse_int(rest[0])
        if len(rest) > 1:
            copy_num = _parse_float(rest[1])
        if len(rest) > 2 and _looks_dna(rest[-1]):
            motif = rest[-1]
        elif len(rest) > 2 and _looks_dna(rest[2]):
            motif = rest[2]
    else:
        # UCSC simpleRepeat: name period copyNum consensusSize ... sequence
        if len(rest) >= 2:
            period = _parse_int(rest[1])
        if len(rest) >= 3:
            copy_num = _parse_float(rest[2])
        if rest and _looks_dna(rest[-1]):
            motif = rest[-1]

    if period is not None and period <= 0:
        period = None
    if motif and period:
        motif = motif[:period]
    if motif:
        motif = motif.upper()

    return chrom, TandemArray(start, end, period, motif, copy_num)


def load_tandem_arrays(bed_path: str) -> dict[str, list[TandemArray]]:
    arrays: dict[str, list[TandemArray]] = defaultdict(list)
    with open_text(bed_path, "rt") as fh:
        for line in fh:
            parsed = parse_tandem_bed_line(line)
            if parsed is None:
                continue
            chrom, arr = parsed
            arrays[chrom].append(arr)
    for chrom in arrays:
        arrays[chrom].sort(key=lambda a: (a.start, a.end))
    return arrays


def overlapping_arrays(
    arrays: list[TandemArray], start0: int, end0: int
) -> list[TandemArray]:
    if end0 <= start0:
        end0 = start0 + 1
    lo, hi = 0, len(arrays)
    while lo < hi:
        mid = (lo + hi) // 2
        if arrays[mid].end <= start0:
            lo = mid + 1
        else:
            hi = mid
    hits = []
    i = lo
    while i < len(arrays) and arrays[i].start < end0:
        a = arrays[i]
        if a.start < end0 and a.end > start0:
            hits.append(a)
        i += 1
    return hits


def infer_svtype(ref: str, alts: list[str], info: dict[str, str]) -> str:
    if "SVTYPE" in info and info["SVTYPE"] not in (".", ""):
        return info["SVTYPE"].split(":")[0].split(",")[0].upper()
    if any("[" in a or "]" in a for a in alts):
        return "BND"
    symbolic = [a.strip("<>").split(":")[0].upper() for a in alts if a.startswith("<")]
    if symbolic:
        return symbolic[0]
    if len(alts) == 1 and len(ref) == 1 and len(alts[0]) == 1 and alts[0] != "*":
        return "SNV"
    if all(len(a) == len(ref) and not a.startswith("<") for a in alts):
        return "MNV"
    svlens = [len(a) - len(ref) for a in alts if not a.startswith("<")]
    if svlens and all(s < 0 for s in svlens):
        return "DEL"
    if svlens and all(s > 0 for s in svlens):
        return "INS"
    return "UNK"


def alt_svlens(
    ref: str, alts: list[str], info: dict[str, str]
) -> list[Optional[int]]:
    raw = info.get("SVLEN", "")
    parsed: list[Optional[int]] = []
    if raw and raw != ".":
        for tok in raw.split(","):
            v = _parse_int(tok)
            parsed.append(v)
    out: list[Optional[int]] = []
    for i, alt in enumerate(alts):
        if i < len(parsed) and parsed[i] is not None:
            out.append(parsed[i])
            continue
        if alt.startswith("<") or "[" in alt or "]" in alt:
            if len(parsed) == 1 and parsed[0] is not None and len(alts) == 1:
                out.append(parsed[0])
            else:
                out.append(None)
            continue
        out.append(len(alt) - len(ref))
    return out


def variant_span(
    pos: int, ref: str, svtype: str, info: dict[str, str], svlens: list[Optional[int]]
) -> tuple[int, int]:
    start0 = pos - 1
    if "END" in info:
        end_v = _parse_int(info["END"])
        if end_v is not None and end_v > start0:
            end0 = end_v
            if svtype == "INS":
                end0 = max(end0, start0 + max(len(ref), 1) + 1)
            return start0, end0
    if svtype == "INS":
        return start0, start0 + max(len(ref), 1) + 1
    abs_lens = [abs(s) for s in svlens if s is not None]
    span = max(abs_lens) if abs_lens else max(len(ref), 1)
    return start0, start0 + max(span, 1)


def indel_sequence(ref: str, alt: str) -> str:
    """Inserted or deleted sequence for a resolved allele."""
    if alt.startswith("<") or "[" in alt or "]" in alt:
        return ""
    if alt.startswith(ref):
        return alt[len(ref) :]
    if ref.startswith(alt):
        return ref[len(alt) :]
    if len(alt) > len(ref):
        return alt
    if len(ref) > len(alt):
        return ref
    return ""


def smallest_period(seq: str) -> Optional[int]:
    n = len(seq)
    if n < 2:
        return None
    for p in range(1, n // 2 + 1):
        if n % p == 0 and seq == seq[:p] * (n // p):
            return p
    return None


def _units_for_period(svlen: Optional[int], seq: str, period: int) -> Optional[int]:
    if period <= 0:
        return None
    if svlen is not None and abs(svlen) % period == 0:
        units = abs(svlen) // period
        if units >= 1:
            return units
    if seq and len(seq) % period == 0:
        units = len(seq) // period
        if units >= 1 and seq == seq[:period] * units:
            return units
    return None


def motif_from_seq(seq: str, period: int, fallback: Optional[str]) -> Optional[str]:
    if fallback:
        return fallback.upper()
    if seq and len(seq) >= period:
        return seq[:period].upper()
    return None


def annotate_record(
    pos: int,
    ref: str,
    alt_field: str,
    info: dict[str, str],
    arrays: list[TandemArray],
) -> Optional[RuAnnotation]:
    alts = [a for a in alt_field.split(",") if a]
    if not alts:
        return None
    svtype = infer_svtype(ref, alts, info)
    if svtype in RU_EXCLUDED or svtype not in RU_ELIGIBLE:
        return None
    svlens = alt_svlens(ref, alts, info)
    start0, end0 = variant_span(pos, ref, svtype, info, svlens)
    hits = overlapping_arrays(arrays, start0, end0)
    if not hits:
        return None

    seqs = [indel_sequence(ref, a) for a in alts]

    def try_period(period: int, motif: Optional[str], cn_ref: Optional[float]) -> Optional[RuAnnotation]:
        rus: list[int] = []
        for svlen, seq in zip(svlens, seqs):
            units = _units_for_period(svlen, seq, period)
            if units is None:
                return None
            rus.append(units)
        use_motif = motif
        if use_motif is None:
            for seq in seqs:
                use_motif = motif_from_seq(seq, period, None)
                if use_motif:
                    break
        return RuAnnotation(period, use_motif, cn_ref, tuple(rus))

    scored: list[tuple[int, int, RuAnnotation]] = []
    for hit in hits:
        if hit.period:
            ann = try_period(hit.period, hit.motif, hit.copy_num)
            if ann is not None:
                overlap = min(hit.end, end0) - max(hit.start, start0)
                scored.append((overlap, -hit.period, ann))
    if scored:
        scored.sort(reverse=True)
        return scored[0][2]

    # 3-column (or period-less) overlap: infer a tiling period from allele sequence.
    inferred: Optional[int] = None
    inferred_motif: Optional[str] = None
    for seq in seqs:
        p = smallest_period(seq)
        if p is None:
            continue
        if inferred is None:
            inferred = p
            inferred_motif = seq[:p].upper()
        elif p != inferred:
            return None
    if inferred is None:
        return None
    cn_ref = next((h.copy_num for h in hits if h.copy_num is not None), None)
    return try_period(inferred, inferred_motif, cn_ref)


def escape_info_string(value: str) -> str:
    return (
        value.replace("%", "%25")
        .replace(";", "%3B")
        .replace("=", "%3D")
        .replace(",", "%2C")
        .replace(" ", "%20")
    )


def format_cn_ref(cn: float) -> str:
    if cn == int(cn):
        return str(int(cn))
    return f"{cn:.6g}"


def merge_info(info: str, ann: Optional[RuAnnotation]) -> str:
    keys = parse_info(info)
    for k in ("RU_TEST", "PERIOD", "MOTIF", "CN_REF", "RU"):
        keys.pop(k, None)
    kept: list[str] = []
    if info and info != ".":
        for item in info.split(";"):
            if not item:
                continue
            key = item.split("=", 1)[0]
            if key in {"RU_TEST", "PERIOD", "MOTIF", "CN_REF", "RU"}:
                continue
            kept.append(item)
    if ann is None:
        return ";".join(kept) if kept else "."
    extra = [
        "RU_TEST",
        f"PERIOD={ann.period}",
        f"RU={','.join(str(u) for u in ann.ru)}",
    ]
    if ann.motif:
        extra.insert(2, f"MOTIF={escape_info_string(ann.motif)}")
    if ann.cn_ref is not None:
        extra.insert(-1, f"CN_REF={format_cn_ref(ann.cn_ref)}")
    return ";".join(kept + extra)


def _info_header_id(line: str) -> Optional[str]:
    if not line.startswith("##INFO=<ID="):
        return None
    rest = line[len("##INFO=<ID=") :]
    return rest.split(",", 1)[0]


def annotate_vcf(in_fh: Iterable[str], out_fh: TextIO, arrays: dict[str, list[TandemArray]]) -> int:
    n_ru = 0
    existing_ids: set[str] = set()
    header_emitted = False
    for line in in_fh:
        if line.startswith("##"):
            iid = _info_header_id(line)
            if iid:
                existing_ids.add(iid)
            out_fh.write(line if line.endswith("\n") else line + "\n")
            continue
        if line.startswith("#"):
            if not header_emitted:
                for hdr in INFO_HEADERS:
                    hid = _info_header_id(hdr)
                    if hid and hid not in existing_ids:
                        out_fh.write(hdr + "\n")
                header_emitted = True
            out_fh.write(line if line.endswith("\n") else line + "\n")
            continue
        if not line.strip():
            continue
        parts = line.rstrip("\n").split("\t")
        if len(parts) < 8:
            out_fh.write(line if line.endswith("\n") else line + "\n")
            continue
        chrom, pos_s, _vid, ref, alt, _qual, _filt, info = parts[:8]
        pos = _parse_int(pos_s)
        ann = None
        if pos is not None:
            ann = annotate_record(pos, ref, alt, parse_info(info), arrays.get(chrom, []))
        if ann is not None:
            n_ru += 1
        parts[7] = merge_info(info, ann)
        out_fh.write("\t".join(parts) + "\n")
    return n_ru


def open_vcf_write(path: str) -> TextIO:
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    if str(path).endswith(".gz"):
        return gzip.open(path, "wt")
    return open(path, "wt")


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", required=True, help="Joint phased VCF (.vcf or .vcf.gz)")
    p.add_argument("--out", required=True, help="Annotated VCF (.vcf or .vcf.gz)")
    p.add_argument(
        "--simple-repeat-bed",
        required=True,
        help="simpleRepeat/TRF BED (3-column or full UCSC/TRF columns)",
    )
    args = p.parse_args()

    arrays = load_tandem_arrays(args.simple_repeat_bed)
    n_arr = sum(len(v) for v in arrays.values())
    print(f"[annotate_ru] loaded {n_arr:,} tandem intervals", file=sys.stderr, flush=True)
    with open_text(args.vcf, "rt") as inf, open_vcf_write(args.out) as out:
        n_ru = annotate_vcf(inf, out, arrays)
    print(f"[annotate_ru] marked {n_ru:,} RU_TEST sites → {args.out}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
