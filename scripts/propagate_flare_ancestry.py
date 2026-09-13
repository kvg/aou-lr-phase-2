#!/usr/bin/env python3
"""Copy FLARE FORMAT/AN1,AN2 onto interstitial sites of a target VCF.

FLARE LAI markers are a sparse subset of the phased callset (sites that overlap
the reference panel and pass MAF/MAC). This fills AN1/AN2 on every target
record from the covering FLARE interval, using the same rules as FELIXla:

* the first FLARE site on a contig covers ``1..current_lai_pos``
* later FLARE sites cover ``(previous_lai_pos, current_lai_pos]``
* the last FLARE state extends through the last target site on that contig

Target ``GT`` and other FORMAT fields are kept. Ancestry tags are matched by
sample name. Streaming is one pass over each sorted VCF (O(sites); one FLARE
record of state).

Example::

    python3 scripts/propagate_flare_ancestry.py \\
      --flare chr22.anc.vcf.gz \\
      --target chr22.phased.vcf.gz \\
      --output chr22.gt_an.vcf \\
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
    "AN1": '##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="Local ancestry of haplotype 1 copied from FLARE">',
    "AN2": '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="Local ancestry of haplotype 2 copied from FLARE">',
    "ANP1": '##FORMAT=<ID=ANP1,Number=1,Type=Float,Description="Ancestry probability of haplotype 1 copied from FLARE">',
    "ANP2": '##FORMAT=<ID=ANP2,Number=1,Type=Float,Description="Ancestry probability of haplotype 2 copied from FLARE">',
}

INFO_FLARE_POS = (
    '##INFO=<ID=FLARE_POS,Number=1,Type=Integer,'
    'Description="POS of the covering FLARE LAI marker (FELIXla interval semantics)">'
)

DEFAULT_TAGS = ("AN1", "AN2")


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
class PropagateStats:
    mode: str = "target_backbone_flare_ancestry"
    region: Optional[str] = None
    tags: list[str] = field(default_factory=lambda: list(DEFAULT_TAGS))
    target_sites: int = 0
    flare_sites: int = 0
    annotated_sites: int = 0
    exact_pos_overlap: int = 0
    sites_before_first_flare: int = 0
    sites_after_last_flare: int = 0
    sites_missing_ancestry: int = 0
    coverage_rate: Optional[float] = None
    n_target_samples: int = 0
    n_flare_samples: int = 0
    n_shared_samples: int = 0
    sample_overlap_rate: Optional[float] = None
    chromosome: str = ""

    def finalize(self) -> None:
        if self.target_sites:
            self.coverage_rate = self.annotated_sites / self.target_sites
        else:
            self.coverage_rate = None
        if self.n_target_samples:
            self.sample_overlap_rate = self.n_shared_samples / self.n_target_samples
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


def _info_has_tag(header_lines: list[str], tag: str) -> bool:
    needle = f"##INFO=<ID={tag},"
    return any(line.startswith(needle) for line in header_lines)


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


def output_format_keys(target_fmt: str, tags: tuple[str, ...]) -> list[str]:
    keys = target_fmt.split(":") if target_fmt else []
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


class FlareCursor:
    """FELIXla covering marker: smallest FLARE POS >= query, else last on contig."""

    def __init__(self, stream: VcfStream):
        self.stream = stream
        self.current: Optional[tuple[str, int, list[str]]] = None
        self.first_pos: dict[str, int] = {}
        self.last_pos: dict[str, int] = {}

    def _skip_to_chrom(self, chrom: str) -> None:
        if self.current is not None and self.current[0] != chrom:
            self.current = None
        while True:
            peek = self.stream.peek()
            if peek is None:
                return
            if peek[0] == chrom:
                return
            self.stream.pop()

    def covering(self, chrom: str, pos: int) -> Optional[tuple[str, int, list[str]]]:
        self._skip_to_chrom(chrom)
        if self.current is None or self.current[0] != chrom:
            rec = self.stream.peek()
            if rec is None or rec[0] != chrom:
                self.current = None
                return None
            self.current = self.stream.pop()
            self.first_pos.setdefault(chrom, self.current[1])
            self.last_pos[chrom] = self.current[1]
            self._skip_same_pos()
        while self.current[1] < pos:
            peek = self.stream.peek()
            if peek is None or peek[0] != chrom:
                return self.current
            self.current = self.stream.pop()
            self.last_pos[chrom] = self.current[1]
            self._skip_same_pos()
        return self.current

    def _skip_same_pos(self) -> None:
        if self.current is None:
            return
        chrom, pos, _ = self.current
        while True:
            peek = self.stream.peek()
            if peek is None or peek[0] != chrom or peek[1] != pos:
                return
            self.stream.pop()


def _extras_for_target(
    flare_cols: Optional[list[str]],
    flare_to_target: list[Optional[int]],
    tags: tuple[str, ...],
) -> list[dict[str, str]]:
    if flare_cols is None:
        return [{} for _ in flare_to_target]
    flare_fmt = flare_cols[8] if len(flare_cols) > 8 else "GT"
    extras: list[dict[str, str]] = []
    for flare_i in flare_to_target:
        if flare_i is None or 9 + flare_i >= len(flare_cols):
            extras.append({})
        else:
            extras.append(extract_tags(flare_fmt, flare_cols[9 + flare_i], tags))
    return extras


def _append_info(info: str, key: str, value: str) -> str:
    extra = f"{key}={value}"
    if not info or info == ".":
        return extra
    return f"{info};{extra}"


def propagate_ancestry(
    flare_vcf: str | Path,
    target_vcf: str | Path,
    output: str | Path | TextIO,
    *,
    tags: tuple[str, ...] = DEFAULT_TAGS,
    region: Optional[Region] = None,
    progress_every: int = 100_000,
    missing_tsv: Optional[Path] = None,
) -> PropagateStats:
    flare = VcfStream(flare_vcf, "flare")
    target = VcfStream(target_vcf, "target")
    stats = PropagateStats(
        tags=list(tags),
        region=None if region is None else f"{region.chrom}:{region.start}-{region.end}",
    )
    try:
        flare_index = {name: i for i, name in enumerate(flare.samples)}
        flare_to_target = [flare_index.get(name) for name in target.samples]
        stats.n_target_samples = len(target.samples)
        stats.n_flare_samples = len(flare.samples)
        stats.n_shared_samples = sum(i is not None for i in flare_to_target)

        if output == "-" or output == "/dev/stdout":
            out_fh: TextIO = sys.stdout
            close_out = False
        elif hasattr(output, "write"):
            out_fh = output  # type: ignore[assignment]
            close_out = False
        else:
            out_fh = open_text(output, "wt")
            close_out = True

        missing_fh: Optional[TextIO] = None
        if missing_tsv is not None:
            missing_tsv.parent.mkdir(parents=True, exist_ok=True)
            missing_fh = missing_tsv.open("w")
            missing_fh.write("chrom\tpos\tref\talt\n")

        cursor = FlareCursor(flare)
        try:
            _write_header(out_fh, target.header, target.samples, flare.header, tags)
            while True:
                rec = target.peek()
                if rec is None:
                    break
                chrom, pos, _ = rec
                if region is not None:
                    if chrom != region.chrom or pos < region.start:
                        target.pop()
                        continue
                    if pos > region.end:
                        break
                cols = target.pop()[2]
                if not stats.chromosome:
                    stats.chromosome = chrom
                stats.target_sites += 1
                covering = cursor.covering(chrom, pos)
                if covering is None:
                    stats.sites_missing_ancestry += 1
                    extras = _extras_for_target(None, flare_to_target, tags)
                    if missing_fh is not None:
                        missing_fh.write(f"{cols[0]}\t{cols[1]}\t{cols[3]}\t{cols[4]}\n")
                    _write_annotated(out_fh, cols, extras, tags, covering_pos=None)
                else:
                    stats.annotated_sites += 1
                    cover_pos = covering[1]
                    if cover_pos == pos:
                        stats.exact_pos_overlap += 1
                    first = cursor.first_pos.get(chrom)
                    last = cursor.last_pos.get(chrom)
                    if first is not None and pos < first:
                        stats.sites_before_first_flare += 1
                    if last is not None and pos > last:
                        stats.sites_after_last_flare += 1
                    extras = _extras_for_target(covering[2], flare_to_target, tags)
                    _write_annotated(out_fh, cols, extras, tags, covering_pos=cover_pos)
                if progress_every and stats.target_sites % progress_every == 0:
                    print(
                        f"[propagate_flare_ancestry] {stats.target_sites:,} target sites "
                        f"({stats.annotated_sites:,} with ancestry)",
                        file=sys.stderr,
                        flush=True,
                    )
        finally:
            if close_out:
                out_fh.close()
            if missing_fh is not None:
                missing_fh.close()
            while flare.peek() is not None:
                flare.pop()
    finally:
        target.close()
        flare.close()
        stats.flare_sites = flare.n_sites
    stats.finalize()
    return stats


def _format_header_for_tag(flare_header: list[str], tag: str) -> str:
    needle = f"##FORMAT=<ID={tag},"
    for line in flare_header:
        if line.startswith(needle):
            return line
    return TAG_HEADERS.get(
        tag,
        f'##FORMAT=<ID={tag},Number=1,Type=String,Description="Copied from FLARE">',
    )


def _write_header(
    out: TextIO,
    target_header: list[str],
    samples: list[str],
    flare_header: list[str],
    tags: tuple[str, ...],
) -> None:
    chrom_line = target_header[-1]
    meta = target_header[:-1]
    for tag in tags:
        if not _format_has_tag(meta, tag):
            meta.append(_format_header_for_tag(flare_header, tag))
    if not _info_has_tag(meta, "FLARE_POS"):
        meta.append(INFO_FLARE_POS)
    meta.append(
        "##propagate_flare_ancestry=FORMAT ancestry tags copied from the covering "
        "FLARE LAI marker (FELIXla interval semantics) by sample name"
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
    covering_pos: Optional[int],
) -> None:
    target_fmt = cols[8] if len(cols) > 8 else "GT"
    new_fmt = output_format_keys(target_fmt, tags)
    n_samp = len(extras_by_sample)
    samples = cols[9 : 9 + n_samp]
    if len(samples) < n_samp:
        samples = samples + ["."] * (n_samp - len(samples))
    new_samples = [
        annotate_sample(target_fmt.split(":") if target_fmt else [], samp, extra, tags)
        for samp, extra in zip(samples, extras_by_sample)
    ]
    info = cols[7] if len(cols) > 7 else "."
    if covering_pos is not None:
        info = _append_info(info, "FLARE_POS", str(covering_pos))
    out_cols = cols[:7] + [info, ":".join(new_fmt)] + new_samples
    out.write("\t".join(out_cols) + "\n")


def _check_thresholds(
    stats: PropagateStats,
    min_coverage_rate: Optional[float],
    min_sample_overlap: Optional[float],
) -> None:
    if min_sample_overlap is not None and stats.sample_overlap_rate is not None:
        if stats.sample_overlap_rate < min_sample_overlap:
            raise SystemExit(
                f"sample overlap too low: {stats.sample_overlap_rate:.4f} "
                f"< {min_sample_overlap:.4f} "
                f"({stats.n_shared_samples}/{stats.n_target_samples} target samples in FLARE)"
            )
    if min_coverage_rate is not None and stats.coverage_rate is not None:
        if stats.coverage_rate < min_coverage_rate:
            raise SystemExit(
                f"FLARE ancestry coverage too low: {stats.coverage_rate:.4f} "
                f"< {min_coverage_rate:.4f} "
                f"({stats.annotated_sites}/{stats.target_sites} target sites)"
            )


def write_stats(stats: PropagateStats, json_path: Optional[Path], tsv_path: Optional[Path]) -> None:
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


def flare_query_region(region: str) -> str:
    """Left-pad a CHR:START-END window so the covering marker before START is kept."""
    spec = region.strip()
    parsed = Region.parse(spec)
    if ":" not in spec or parsed.end == sys.maxsize:
        return parsed.chrom
    return f"{parsed.chrom}:1-{parsed.end}"


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--flare", required=True, help="FLARE ancestry VCF (.vcf / .vcf.gz)")
    p.add_argument(
        "--target",
        required=True,
        help="Target backbone VCF (phased SNV/indel, SV, or joint). GT is kept.",
    )
    p.add_argument("--output", required=True, help="Output VCF path, or - for stdout")
    p.add_argument("--stats-json", type=Path)
    p.add_argument("--stats-tsv", type=Path)
    p.add_argument(
        "--missing-tsv",
        type=Path,
        help="Target sites with no FLARE marker on the contig",
    )
    p.add_argument(
        "--tags",
        default="AN1,AN2",
        help="Comma-separated FLARE FORMAT tags to copy (default: AN1,AN2)",
    )
    p.add_argument("--region", help="Optional CHR or CHR:START-END filter on the target")
    p.add_argument(
        "--min-coverage-rate",
        type=float,
        default=None,
        help="Fail if annotated target sites / target sites is below this fraction",
    )
    p.add_argument(
        "--min-sample-overlap",
        type=float,
        default=None,
        help="Fail if shared samples / target samples is below this fraction",
    )
    p.add_argument("--progress-every", type=int, default=100_000)
    args = p.parse_args(argv)

    tags = tuple(t.strip() for t in args.tags.split(",") if t.strip())
    if not tags:
        raise SystemExit("--tags must list at least one FORMAT field")
    region = Region.parse(args.region) if args.region else None
    stats = propagate_ancestry(
        args.flare,
        args.target,
        args.output,
        tags=tags,
        region=region,
        progress_every=args.progress_every,
        missing_tsv=args.missing_tsv,
    )
    write_stats(stats, args.stats_json, args.stats_tsv)
    rate = "NA" if stats.coverage_rate is None else f"{stats.coverage_rate * 100:.2f}%"
    print(
        f"OK: {stats.annotated_sites}/{stats.target_sites} target sites received FLARE ancestry "
        f"({rate}; exact_pos={stats.exact_pos_overlap} missing={stats.sites_missing_ancestry} "
        f"before_first={stats.sites_before_first_flare} after_last={stats.sites_after_last_flare}); "
        f"{stats.n_shared_samples}/{stats.n_target_samples} samples shared",
        file=sys.stderr,
    )
    _check_thresholds(stats, args.min_coverage_rate, args.min_sample_overlap)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
