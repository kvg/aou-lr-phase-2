#!/usr/bin/env python3
"""Copy FORMAT/GQ, FORMAT/DP, and FORMAT/RNC from a DeepVariant VCF onto FLARE sites.

FLARE is the backbone: every FLARE record is emitted with its original GT/AN*
fields. Matching DeepVariant (GLnexus) sites contribute GQ/DP/RNC by sample name.
Streaming is one pass over each sorted VCF (O(sites); one position buffered).

GLnexus ``RNC`` (Reason for No Call, two characters) distinguishes placeholder
``0/0`` (``..``) from incomplete gVCF coverage (``I``).

Example::

    python3 scripts/annotate_flare_gq_dp.py \\
      --flare chr22.anc.vcf.gz \\
      --deepvariant chr22.vcf.gz \\
      --output chr22.flare.gq_dp.vcf \\
      --stats-json stats.json
"""

from __future__ import annotations

import argparse
import gzip
import json
import sys
from dataclasses import asdict, dataclass, field
from pathlib import Path
from typing import Optional, TextIO


TAG_HEADERS = {
    "GQ": '##FORMAT=<ID=GQ,Number=1,Type=Integer,Description="Genotype quality copied from DeepVariant/GLnexus">',
    "DP": '##FORMAT=<ID=DP,Number=1,Type=Integer,Description="Read depth copied from DeepVariant/GLnexus">',
    "RNC": '##FORMAT=<ID=RNC,Number=1,Type=String,Description="GLnexus Reason for No Call copied from DeepVariant/GLnexus (. or .. = called; I = incomplete gVCF)">',
}

DEFAULT_TAGS = ("GQ", "DP", "RNC")


@dataclass(frozen=True)
class Region:
    chrom: str
    start: int
    end: int

    @classmethod
    def parse(cls, spec: str) -> "Region":
        spec = spec.strip()
        if not spec:
            raise ValueError("empty region")
        if ":" not in spec:
            return cls(chrom=spec, start=1, end=sys.maxsize)
        chrom, rest = spec.split(":", 1)
        if not chrom:
            raise ValueError(f"empty chromosome in region {spec!r}")
        if "-" not in rest:
            raise ValueError(f"expected CHR:START-END, got {spec!r}")
        start_s, end_s = rest.split("-", 1)
        start = int(start_s)
        end = sys.maxsize if end_s == "" else int(end_s)
        if start < 1 or end < start:
            raise ValueError(f"start must be >= 1 and end >= start: {spec!r}")
        return cls(chrom=chrom, start=start, end=end)

    def contains(self, chrom: str, pos: int) -> bool:
        return chrom == self.chrom and self.start <= pos <= self.end


@dataclass
class AnnotateStats:
    mode: str = "flare_backbone_gq_dp"
    region: Optional[str] = None
    tags: list[str] = field(default_factory=lambda: list(DEFAULT_TAGS))
    flare_sites: int = 0
    deepvariant_sites: int = 0
    matched_sites: int = 0
    matched_exact: int = 0
    matched_multiallelic: int = 0
    unmatched_flare_sites: int = 0
    match_rate: Optional[float] = None
    n_flare_samples: int = 0
    n_deepvariant_samples: int = 0
    n_shared_samples: int = 0
    sample_overlap_rate: Optional[float] = None
    chromosome: str = ""

    def finalize(self) -> None:
        if self.flare_sites:
            self.match_rate = self.matched_sites / self.flare_sites
        else:
            self.match_rate = None
        if self.n_flare_samples:
            self.sample_overlap_rate = self.n_shared_samples / self.n_flare_samples
        else:
            self.sample_overlap_rate = None


def open_text(path: str | Path, mode: str = "rt"):
    path = Path(path)
    if str(path).endswith(".gz") or str(path).endswith(".bgz"):
        return gzip.open(path, mode)
    return open(path, mode)


def _format_has_tag(header_lines: list[str], tag: str) -> bool:
    needle = f"##FORMAT=<ID={tag},"
    return any(line.startswith(needle) for line in header_lines)


def match_dv_record(
    ref: str, alt: str, dv_recs: list[list[str]]
) -> tuple[Optional[list[str]], str]:
    """Match a FLARE allele to a DeepVariant record at the same POS.

    ``exact``: same REF and ALT strings.
    ``multiallelic``: same REF, and at least one FLARE ALT is in the DV ALT list
    (GLnexus often keeps A,B while FLARE is split, or the reverse).
    """
    exact: dict[tuple[str, str], list[str]] = {}
    same_ref: list[list[str]] = []
    for rec in dv_recs:
        exact[(rec[3], rec[4])] = rec
        if rec[3] == ref:
            same_ref.append(rec)
    hit = exact.get((ref, alt))
    if hit is not None:
        return hit, "exact"
    flare_alts = {a for a in alt.split(",") if a and a != "."}
    if not flare_alts:
        return None, ""
    for rec in same_ref:
        dv_alts = {a for a in rec[4].split(",") if a and a != "."}
        if flare_alts & dv_alts:
            return rec, "multiallelic"
    return None, ""


def extract_tags(fmt: str, sample: str, tags: tuple[str, ...]) -> dict[str, str]:
    keys = fmt.split(":")
    vals = sample.split(":")
    idx = {k: i for i, k in enumerate(keys)}
    out: dict[str, str] = {}
    for tag in tags:
        i = idx.get(tag)
        if i is None or i >= len(vals) or vals[i] == "":
            out[tag] = "."
        else:
            out[tag] = vals[i]
    return out


def annotate_sample(
    fmt_keys: list[str],
    sample: str,
    extra: dict[str, str],
    tags: tuple[str, ...],
) -> str:
    vals = sample.split(":")
    if len(vals) < len(fmt_keys):
        vals = vals + ["."] * (len(fmt_keys) - len(vals))
    else:
        vals = vals[: len(fmt_keys)]
    key_idx = {k: i for i, k in enumerate(fmt_keys)}
    out_vals = list(vals)
    for tag in tags:
        value = extra.get(tag, ".")
        if tag in key_idx:
            out_vals[key_idx[tag]] = value
        else:
            out_vals.append(value)
    return ":".join(out_vals)


def output_format_keys(flare_fmt: str, tags: tuple[str, ...]) -> list[str]:
    keys = flare_fmt.split(":") if flare_fmt else []
    for tag in tags:
        if tag not in keys:
            keys.append(tag)
    return keys


class VcfStream:
    """Sorted VCF record stream. Buffers at most one record plus the current POS."""

    def __init__(self, path: str | Path, label: str):
        self.path = Path(path)
        self.label = label
        self._fh: TextIO = open_text(self.path)
        self.header: list[str] = []
        self.samples: list[str] = []
        self._peek: Optional[tuple[str, int, list[str]]] = None
        self._prev: Optional[tuple[str, int]] = None
        self.n_sites = 0
        self._read_header()

    def _read_header(self) -> None:
        for line in self._fh:
            if line.startswith("##"):
                self.header.append(line.rstrip("\n"))
                continue
            if line.startswith("#CHROM"):
                cols = line.rstrip("\n").split("\t")
                if len(cols) < 9:
                    raise ValueError(f"{self.label}: malformed #CHROM line")
                self.samples = cols[9:]
                self.header.append(line.rstrip("\n"))
                return
            raise ValueError(f"{self.label}: missing #CHROM header")
        raise ValueError(f"{self.label}: empty or header-only VCF")

    def close(self) -> None:
        self._fh.close()

    def _read_one(self) -> Optional[tuple[str, int, list[str]]]:
        while True:
            line = self._fh.readline()
            if not line:
                return None
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 8:
                raise ValueError(f"{self.label}: short VCF line: {line[:80]!r}")
            chrom = cols[0]
            try:
                pos = int(cols[1])
            except ValueError as exc:
                raise ValueError(f"{self.label}: bad POS {cols[1]!r}") from exc
            key = (chrom, pos)
            if self._prev is not None and key < self._prev:
                raise ValueError(
                    f"{self.label}: unsorted VCF: {key} follows {self._prev}"
                )
            self._prev = key
            self.n_sites += 1
            return chrom, pos, cols

    def peek(self) -> Optional[tuple[str, int, list[str]]]:
        if self._peek is None:
            self._peek = self._read_one()
        return self._peek

    def pop(self) -> Optional[tuple[str, int, list[str]]]:
        rec = self.peek()
        self._peek = None
        return rec

    def take_position(self, chrom: str, pos: int) -> list[list[str]]:
        recs: list[list[str]] = []
        while True:
            cur = self.peek()
            if cur is None:
                break
            c, p, cols = cur
            if c != chrom or p != pos:
                break
            recs.append(cols)
            self.pop()
        return recs


def _skip_before(stream: VcfStream, chrom: str, pos: int) -> None:
    while True:
        cur = stream.peek()
        if cur is None:
            return
        c, p, _ = cur
        if (c, p) >= (chrom, pos):
            return
        stream.pop()


def annotate_flare(
    flare_vcf: str | Path,
    deepvariant_vcf: str | Path,
    output: str | Path | TextIO,
    *,
    tags: tuple[str, ...] = DEFAULT_TAGS,
    region: Optional[Region] = None,
    progress_every: int = 100_000,
    unmatched_tsv: Optional[Path] = None,
) -> AnnotateStats:
    flare = VcfStream(flare_vcf, "flare")
    dv = VcfStream(deepvariant_vcf, "deepvariant")
    stats = AnnotateStats(tags=list(tags), region=None if region is None else f"{region.chrom}:{region.start}-{region.end}")
    try:
        dv_index = {name: i for i, name in enumerate(dv.samples)}
        flare_to_dv = [dv_index.get(name) for name in flare.samples]
        stats.n_flare_samples = len(flare.samples)
        stats.n_deepvariant_samples = len(dv.samples)
        stats.n_shared_samples = sum(i is not None for i in flare_to_dv)

        if output == "-" or output == "/dev/stdout":
            out_fh: TextIO = sys.stdout
            close_out = False
        elif hasattr(output, "write"):
            out_fh = output  # type: ignore[assignment]
            close_out = False
        else:
            out_fh = open_text(output, "wt")
            close_out = True

        unmatched_fh: Optional[TextIO] = None
        if unmatched_tsv is not None:
            unmatched_tsv.parent.mkdir(parents=True, exist_ok=True)
            unmatched_fh = unmatched_tsv.open("w")
            unmatched_fh.write("chrom\tpos\tref\talt\n")

        try:
            _write_header(out_fh, flare.header, flare.samples, tags)
            while True:
                rec = flare.peek()
                if rec is None:
                    break
                chrom, pos, _ = rec
                if region is not None:
                    if chrom != region.chrom or pos < region.start:
                        flare.pop()
                        continue
                    if pos > region.end:
                        break
                if not stats.chromosome:
                    stats.chromosome = chrom
                _skip_before(dv, chrom, pos)
                flare_recs = flare.take_position(chrom, pos)
                dv_at = dv.peek()
                dv_recs: list[list[str]] = []
                if dv_at is not None and dv_at[0] == chrom and dv_at[1] == pos:
                    dv_recs = dv.take_position(chrom, pos)

                for cols in flare_recs:
                    stats.flare_sites += 1
                    match, how = match_dv_record(cols[3], cols[4], dv_recs)
                    if match is None:
                        stats.unmatched_flare_sites += 1
                        extras_by_sample = [{} for _ in flare.samples]
                        if unmatched_fh is not None:
                            unmatched_fh.write(
                                f"{cols[0]}\t{cols[1]}\t{cols[3]}\t{cols[4]}\n"
                            )
                    else:
                        stats.matched_sites += 1
                        if how == "exact":
                            stats.matched_exact += 1
                        else:
                            stats.matched_multiallelic += 1
                        dv_fmt = match[8] if len(match) > 8 else "GT"
                        extras_by_sample = []
                        for dv_i in flare_to_dv:
                            if dv_i is None or 9 + dv_i >= len(match):
                                extras_by_sample.append({})
                            else:
                                extras_by_sample.append(
                                    extract_tags(dv_fmt, match[9 + dv_i], tags)
                                )
                    _write_annotated(out_fh, cols, extras_by_sample, tags)
                    if progress_every and stats.flare_sites % progress_every == 0:
                        print(
                            f"[annotate_flare_gq_dp] {stats.flare_sites:,} FLARE sites "
                            f"({stats.matched_sites:,} matched)",
                            file=sys.stderr,
                            flush=True,
                        )
        finally:
            if close_out:
                out_fh.close()
            if unmatched_fh is not None:
                unmatched_fh.close()
    finally:
        flare.close()
        dv.close()

        stats.deepvariant_sites = dv.n_sites
    stats.finalize()
    return stats


def _write_header(
    out: TextIO,
    header: list[str],
    samples: list[str],
    tags: tuple[str, ...],
) -> None:
    chrom_line = header[-1]
    meta = header[:-1]
    for tag in tags:
        if not _format_has_tag(meta, tag):
            meta.append(TAG_HEADERS.get(tag, f'##FORMAT=<ID={tag},Number=1,Type=String,Description="Copied from DeepVariant">'))
    meta.append(
        "##annotate_flare_gq_dp=FORMAT/GQ,FORMAT/DP,FORMAT/RNC copied from DeepVariant/GLnexus by sample name"
    )
    for line in meta:
        out.write(line + "\n")
    cols = chrom_line.split("\t")
    cols = cols[:9] + samples
    out.write("\t".join(cols) + "\n")


def _write_annotated(
    out: TextIO,
    cols: list[str],
    extras_by_sample: list[dict[str, str]],
    tags: tuple[str, ...],
) -> None:
    flare_fmt = cols[8] if len(cols) > 8 else "GT"
    new_fmt = output_format_keys(flare_fmt, tags)
    n_samp = len(extras_by_sample)
    samples = cols[9 : 9 + n_samp]
    if len(samples) < n_samp:
        samples = samples + ["."] * (n_samp - len(samples))
    new_samples = [
        annotate_sample(flare_fmt.split(":") if flare_fmt else [], samp, extra, tags)
        for samp, extra in zip(samples, extras_by_sample)
    ]
    out_cols = cols[:8] + [":".join(new_fmt)] + new_samples
    out.write("\t".join(out_cols) + "\n")


def _check_thresholds(
    stats: AnnotateStats,
    min_match_rate: Optional[float],
    min_sample_overlap: Optional[float],
) -> None:
    if min_sample_overlap is not None and stats.sample_overlap_rate is not None:
        if stats.sample_overlap_rate < min_sample_overlap:
            raise SystemExit(
                f"sample overlap too low: {stats.sample_overlap_rate:.4f} "
                f"< {min_sample_overlap:.4f} "
                f"({stats.n_shared_samples}/{stats.n_flare_samples} FLARE samples in DeepVariant)"
            )
    if min_match_rate is not None and stats.match_rate is not None:
        if stats.match_rate < min_match_rate:
            raise SystemExit(
                f"FLARE–DeepVariant site match rate too low: {stats.match_rate:.4f} "
                f"< {min_match_rate:.4f} "
                f"({stats.matched_sites}/{stats.flare_sites} sites)"
            )


def write_stats(stats: AnnotateStats, json_path: Optional[Path], tsv_path: Optional[Path]) -> None:
    payload = asdict(stats)
    if json_path is not None:
        json_path.parent.mkdir(parents=True, exist_ok=True)
        json_path.write_text(json.dumps(payload, indent=2) + "\n")
    if tsv_path is not None:
        tsv_path.parent.mkdir(parents=True, exist_ok=True)
        with tsv_path.open("w") as fh:
            for key, value in payload.items():
                if isinstance(value, list):
                    value = ",".join(str(v) for v in value)
                fh.write(f"{key}\t{value}\n")


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flare", required=True, help="FLARE ancestry VCF (.vcf / .vcf.gz)")
    p.add_argument(
        "--deepvariant",
        required=True,
        help="DeepVariant/GLnexus joint VCF with FORMAT/GQ, FORMAT/DP, and FORMAT/RNC",
    )
    p.add_argument("--output", required=True, help="Output VCF path, or - for stdout")
    p.add_argument("--stats-json", type=Path)
    p.add_argument("--stats-tsv", type=Path)
    p.add_argument("--unmatched-tsv", type=Path, help="FLARE sites with no DeepVariant GQ/DP/RNC match")
    p.add_argument(
        "--tags",
        default="GQ,DP,RNC",
        help="Comma-separated FORMAT tags to copy (default: GQ,DP,RNC)",
    )
    p.add_argument("--region", help="Optional CHR or CHR:START-END filter")
    p.add_argument(
        "--min-match-rate",
        type=float,
        default=None,
        help="Fail if matched FLARE sites / FLARE sites is below this fraction",
    )
    p.add_argument(
        "--min-sample-overlap",
        type=float,
        default=None,
        help="Fail if shared samples / FLARE samples is below this fraction",
    )
    p.add_argument("--progress-every", type=int, default=100_000)
    args = p.parse_args(argv)

    tags = tuple(t.strip() for t in args.tags.split(",") if t.strip())
    if not tags:
        raise SystemExit("--tags must list at least one FORMAT field")
    region = Region.parse(args.region) if args.region else None
    stats = annotate_flare(
        args.flare,
        args.deepvariant,
        args.output,
        tags=tags,
        region=region,
        progress_every=args.progress_every,
        unmatched_tsv=args.unmatched_tsv,
    )
    write_stats(stats, args.stats_json, args.stats_tsv)
    rate = "NA" if stats.match_rate is None else f"{stats.match_rate * 100:.2f}%"
    print(
        f"OK: {stats.matched_sites}/{stats.flare_sites} FLARE sites matched DeepVariant "
        f"({rate}; exact={stats.matched_exact} multiallelic={stats.matched_multiallelic} "
        f"unmatched={stats.unmatched_flare_sites}); "
        f"{stats.n_shared_samples}/{stats.n_flare_samples} samples shared",
        file=sys.stderr,
    )
    _check_thresholds(stats, args.min_match_rate, args.min_sample_overlap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
