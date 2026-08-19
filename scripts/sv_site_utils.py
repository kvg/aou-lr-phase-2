#!/usr/bin/env python3
"""Shared helpers for SV annotation site tables."""

from __future__ import annotations

import gzip
import math
import os
import shutil
import subprocess
import sys
from contextlib import contextmanager
from pathlib import Path
from typing import Iterator, Optional, TextIO


SITE_COLUMNS = [
    "chrom",
    "pos",
    "end",
    "id",
    "ref",
    "alt",
    "svtype",
    "svlen",
    "filter",
    "source_vcf",
    "phase",
    "ac",
    "an",
    "af",
    "n_carriers",
    "freq_class",
    "region_class",
    "hit_rmsk",
    "hit_simpleRepeat",
    "hit_genomicSuperDups",
    "cadd_sv_phred",
    "cadd_sv_bin",
    "size_bin_ge20",
    "size_bin_ge50",
    "suppressed_main_duplicate",
]


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def file_magic(path: str | Path, n: int = 4) -> bytes:
    with open(path, "rb") as fh:
        return fh.read(n)


def is_binary_variant_file(path: str | Path) -> bool:
    """True for BCF or gzip/bgzip (including .vcf.gz / compressed .bcf)."""
    magic = file_magic(path)
    if magic[:2] == b"\x1f\x8b":
        return True
    if magic[:3] == b"BCF":
        return True
    name = str(path).lower()
    return name.endswith((".bcf", ".bcf.gz", ".vcf.gz", ".vcf.bgz"))


def bcftools_bin() -> Optional[str]:
    return shutil.which("bcftools")


def bcftools_thread_args() -> list[str]:
    n = os.cpu_count() or 1
    return ["--threads", str(max(1, min(int(n), 8)))]


@contextmanager
def bcftools_stdout(args: list[str], path: str) -> Iterator[TextIO]:
    exe = bcftools_bin()
    if not exe:
        raise RuntimeError(
            f"bcftools is required ({path}). The AnnotateSvCallset image includes bcftools."
        )
    proc = subprocess.Popen(
        [exe, *args, path],
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        text=True,
        bufsize=1 << 20,
    )
    assert proc.stdout is not None
    try:
        yield proc.stdout
    finally:
        proc.stdout.close()
        proc.wait()
        if proc.returncode not in (0, None):
            raise RuntimeError(
                f"bcftools {' '.join(args)} failed for {path} (exit {proc.returncode})"
            )


@contextmanager
def open_vcf_text(path: str | Path) -> Iterator[TextIO]:
    """Yield a text VCF stream (full records, including genotypes).

    Prefer ``iter_site_records`` / ``iter_gt_records`` for partition processing;
    those avoid expanding megabase ALT alleles and FORMAT fields.
    """
    path_s = str(path)
    exe = bcftools_bin()
    binary = is_binary_variant_file(path_s)
    if binary or exe:
        if not exe:
            raise RuntimeError(
                f"bcftools is required to read compressed/BCF input ({path_s}). "
                "The AnnotateSvCallset image includes bcftools."
            )
        with bcftools_stdout(["view", "--no-version", *bcftools_thread_args()], path_s) as fh:
            yield fh
        return
    with open_text(path_s, "rt") as fh:
        yield fh


def list_vcf_samples(path: str | Path) -> list[str]:
    path_s = str(path)
    exe = bcftools_bin()
    if exe:
        out = subprocess.check_output([exe, "query", "-l", path_s], text=True)
        return [ln for ln in out.splitlines() if ln]
    with open_vcf_text(path_s) as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                return line.rstrip("\n").split("\t")[9:]
    return []


def iter_site_records(path: str | Path) -> Iterator[tuple[str, str, str, str, str, str, str]]:
    """Yield (chrom, pos, id, ref, alt, filter, info) with no sample columns.

    ``bcftools query`` omits ALT so ultralong insertions are not pulled into
    Python. Missing ALT is ``.``; callers should take SVTYPE/SVLEN from INFO.
    """
    path_s = str(path)
    if bcftools_bin():
        fmt = r"%CHROM\t%POS\t%ID\t%REF\t%FILTER\t%INFO\n"
        with bcftools_stdout(
            ["query", "-f", fmt],
            path_s,
        ) as fh:
            for line in fh:
                chrom, pos, vid, ref, filt, info = line.rstrip("\n").split("\t", 5)
                yield chrom, pos, vid, ref, ".", filt, info
        return
    with open_vcf_text(path_s) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t", 8)
            if len(parts) < 8:
                continue
            chrom, pos, vid, ref, alt, _qual, filt, info = parts[:8]
            if len(alt) > 64 and not (alt.startswith("<") and alt.endswith(">")):
                alt = "."
            yield chrom, pos, vid, ref, alt, filt, info


def iter_gt_records(path: str | Path) -> Iterator[tuple[str, int, str, list[str]]]:
    """Yield (chrom, pos, id, gt_values) — GT only, no REF/ALT/FORMAT extras."""
    path_s = str(path)
    if bcftools_bin():
        fmt = r"%CHROM\t%POS\t%ID[\t%GT]\n"
        with bcftools_stdout(
            ["query", "-f", fmt],
            path_s,
        ) as fh:
            for line in fh:
                parts = line.rstrip("\n").split("\t", 3)
                if len(parts) < 3:
                    continue
                chrom, pos, vid = parts[0], parts[1], parts[2]
                rest = parts[3] if len(parts) > 3 else ""
                gts = rest.split("\t") if rest else []
                yield chrom, int(pos), vid, gts
        return
    with open_vcf_text(path_s) as fh:
        for line in fh:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            parts = line.rstrip("\n").split("\t", 9)
            if len(parts) < 9:
                continue
            chrom, pos, vid, _ref, _alt, _qual, _filt, _info, fmt = parts[:9]
            rest = parts[9] if len(parts) > 9 else ""
            fmt_fields = fmt.split(":")
            try:
                gt_i = fmt_fields.index("GT")
            except ValueError:
                continue
            gts: list[str] = []
            for cell in rest.split("\t"):
                bits = cell.split(":")
                gts.append(bits[gt_i] if gt_i < len(bits) else ".")
            yield chrom, int(pos), vid, gts


def parse_info(info: str) -> dict[str, str]:
    out: dict[str, str] = {}
    if not info or info == ".":
        return out
    for item in info.split(";"):
        if not item:
            continue
        if "=" in item:
            k, v = item.split("=", 1)
            out[k] = v
        else:
            out[item] = "True"
    return out


def abs_svlen(svlen: Optional[int]) -> Optional[int]:
    if svlen is None:
        return None
    return abs(int(svlen))


def infer_svlen(ref: str, alt: str, info: dict[str, str]) -> Optional[int]:
    if "SVLEN" in info and info["SVLEN"] not in (".", ""):
        try:
            return int(info["SVLEN"].split(",")[0])
        except ValueError:
            pass
    if alt.startswith("<") or ref.startswith("<"):
        if "END" in info:
            try:
                end = int(info["END"])
                # symbolic; length unknown without POS in this helper
                return None
            except ValueError:
                return None
        return None
    # sequence-resolved indel-like alleles
    if "," in alt:
        alt = alt.split(",", 1)[0]
    return len(alt) - len(ref)


def infer_svtype(alt: str, info: dict[str, str], source_vcf: str) -> str:
    if source_vcf == "bnd":
        return info.get("SVTYPE", "BND")
    if "SVTYPE" in info and info["SVTYPE"] not in (".", ""):
        return info["SVTYPE"].split(":")[0]
    if alt.startswith("<") and alt.endswith(">"):
        return alt.strip("<>").split(":")[0]
    if "[" in alt or "]" in alt:
        return "BND"
    svlen = infer_svlen("N", alt, info)
    if svlen is None:
        return "UNK"
    if svlen < 0:
        return "DEL"
    if svlen > 0:
        return "INS"
    return "UNK"


def freq_class(ac: int, an: int, n_carriers: int, n_samples: int) -> str:
    if an <= 0:
        af = 0.0
    else:
        af = ac / an
    if n_samples > 0 and n_carriers >= n_samples:
        return "shared"
    if ac <= 1:
        return "singleton"
    if af >= 0.5:
        return "major"
    return "polymorphic"


def cadd_sv_bin(phred: Optional[float]) -> str:
    if phred is None or (isinstance(phred, float) and math.isnan(phred)):
        return "unscored"
    if phred < 10:
        return "low"
    if phred < 20:
        return "mid"
    return "high"


def size_flags(svlen: Optional[int], svtype: str) -> tuple[bool, bool]:
    if svtype == "BND" or svlen is None:
        return False, False
    a = abs(int(svlen))
    return a >= 20, a >= 50


def write_site_header(fh: TextIO) -> None:
    fh.write("\t".join(SITE_COLUMNS) + "\n")


def site_row(values: dict) -> str:
    out = []
    for col in SITE_COLUMNS:
        v = values.get(col, "")
        if v is None:
            v = ""
        elif isinstance(v, bool):
            v = "true" if v else "false"
        out.append(str(v))
    return "\t".join(out) + "\n"


def parse_site_row(header: list[str], line: str) -> dict:
    parts = line.rstrip("\n").split("\t")
    row = dict(zip(header, parts))
    for k in ("pos", "end", "svlen", "ac", "an", "n_carriers"):
        if k in row and row[k] not in ("", "."):
            try:
                row[k] = int(float(row[k]))
            except ValueError:
                row[k] = None
        elif k in row:
            row[k] = None
    if "af" in row and row["af"] not in ("", "."):
        try:
            row["af"] = float(row["af"])
        except ValueError:
            row["af"] = float("nan")
    if "cadd_sv_phred" in row and row["cadd_sv_phred"] not in ("", "."):
        try:
            row["cadd_sv_phred"] = float(row["cadd_sv_phred"])
        except ValueError:
            row["cadd_sv_phred"] = None
    else:
        row["cadd_sv_phred"] = None
    for k in (
        "hit_rmsk",
        "hit_simpleRepeat",
        "hit_genomicSuperDups",
        "size_bin_ge20",
        "size_bin_ge50",
        "suppressed_main_duplicate",
    ):
        if k in row:
            row[k] = str(row[k]).lower() in ("1", "true", "t", "yes")
    return row


def iter_sites(path: str | Path) -> Iterator[dict]:
    """Yield site-table rows without loading the full table into memory."""
    with open_text(path, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")
        for line in fh:
            if not line.strip():
                continue
            yield parse_site_row(header, line)


def read_sites(path: str | Path) -> list[dict]:
    return list(iter_sites(path))


def ancestry_sort_key(label: str) -> tuple[int, int, str]:
    """Return (afr_group, within_group_rank, label). Non-African first."""
    lab = (label or "oth").strip().lower()
    non_afr = ["eas", "amr", "eur", "sas", "oth", "other", "unknown", "nan", ""]
    if lab in ("afr", "african"):
        return (1, 0, lab)
    if lab in non_afr:
        return (0, non_afr.index(lab), lab)
    return (0, 50, lab)


def order_samples_by_ancestry(sample_ancestry: dict[str, str]) -> list[str]:
    return sorted(
        sample_ancestry.keys(),
        key=lambda s: (*ancestry_sort_key(sample_ancestry.get(s, "oth")), s),
    )
