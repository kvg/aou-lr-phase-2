#!/usr/bin/env python3
"""Detect FLARE ancestry switches, tract lengths, and optional GQ / DP QC.

Streams a FLARE VCF/BCF (FORMAT ``AN1``, ``AN2`` required). ``GQ``, ``DP``,
``GT``, and ``RNC`` are used when present (``AnnotateFlareGqDp`` output).
Raw ``FlareByPopulation`` ``*.anc.vcf.gz`` files work for switch rates and
tract lengths. Inputs must be sorted by chromosome position.

GQ/DP are stratified by genotype class (hom_ref / het / hom_alt) so the
GLnexus homozygous-reference sentinel (typically ``0/0`` with DP=0, GQ≤1) is
not mixed with low-quality hets.

Most FLARE switches land on ``0/0``. For those events the genotype at the
switch site is uninformative; the scan also records GQ/DP of the previous het
and the nearest het (previous or next in bp) on the same sample.

Tract lengths are the ancestry stretches between switches, plus window-edge
segments. In a 1 Mb slice most tracts are censored; **switches per hap per Mb**
is the metric to compare against T (≈ T/100 at 1 cM/Mb).

Example::

    python3 scripts/flare_switch_qc.py \\
      --vcf annotated.chr22.bcf \\
      --samples analysis_samples.txt \\
      --region chr22:26897597-27897597 \\
      --out-dir flare_switch_qc_chr22 \\
      --gq-threshold 20 \\
      --dp-threshold 5

For association-grid scoring (rates / flicker only), prefer::

    python3 scripts/flare_switch_qc.py \\
      --vcf run.anc.vcf.gz --region chr22:start-end \\
      --out-dir out --an-only --summary-only --jobs 8
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import re
import shutil
import subprocess
import sys
from collections import Counter, defaultdict
from concurrent.futures import ProcessPoolExecutor
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional, Sequence, TextIO


ANCESTRY = {0: "eas", 1: "amr", 2: "eur", 3: "afr", 4: "sas"}
GT_CLASSES = ("hom_ref", "het", "hom_alt", "missing", "other")
# Admixture-tract expectation at 1 cM/Mb: mean length Mb = 100/T, rate = T/100.
DEFAULT_REC_CM_PER_MB = 1.0


@dataclass
class HapState:
    anc: Optional[int] = None
    pos: Optional[int] = None
    gq: Optional[int] = None
    dp: Optional[int] = None


@dataclass
class SwitchEvent:
    sample: str
    hap: int  # 1 or 2
    chrom: str
    pos: int
    ref: str
    alt: str
    anc_from: int
    anc_to: int
    gq: Optional[int]
    dp: Optional[int]
    rnc: str
    prev_pos: int
    bp_gap: int
    is_flicker: bool = False
    flicker_span_bp: Optional[int] = None
    gt: str = "."
    gt_class: str = "missing"
    is_glnexus_homref_sentinel: bool = False
    prev_het_pos: Optional[int] = None
    prev_het_gq: Optional[int] = None
    prev_het_dp: Optional[int] = None
    prev_het_gap: Optional[int] = None
    next_het_pos: Optional[int] = None
    next_het_gq: Optional[int] = None
    next_het_dp: Optional[int] = None
    next_het_gap: Optional[int] = None
    nearest_het_pos: Optional[int] = None
    nearest_het_gq: Optional[int] = None
    nearest_het_dp: Optional[int] = None
    nearest_het_gap: Optional[int] = None
    nearest_het_side: str = "none"
    prev_het_rnc: str = "."
    next_het_rnc: str = "."
    nearest_het_rnc: str = "."


@dataclass
class AncTract:
    """Ancestry stretch on one hap, possibly censored at the scan window."""

    sample: str
    hap: int
    chrom: str
    start: int
    end: int
    anc: Optional[int]
    left_censored: bool
    right_censored: bool
    n_switches_on_hap: int
    is_flicker_interval: bool = False

    @property
    def length_bp(self) -> int:
        return max(0, self.end - self.start)


@dataclass
class SampleTracker:
    hap1: HapState = field(default_factory=HapState)
    hap2: HapState = field(default_factory=HapState)
    switches: list[SwitchEvent] = field(default_factory=list)
    last_chrom: Optional[str] = None
    last_het_pos: Optional[int] = None
    last_het_gq: Optional[int] = None
    last_het_dp: Optional[int] = None
    last_het_rnc: str = "."
    pending_next_het: list[SwitchEvent] = field(default_factory=list)


@dataclass
class GtBin:
    """Per-genotype-class quality accumulator."""

    n: int = 0
    low_gq: int = 0
    low_dp: int = 0
    sentinel: int = 0
    rnc_i: int = 0
    gq: list[int] = field(default_factory=list)
    dp: list[int] = field(default_factory=list)


@dataclass
class QcAccum:
    """On-the-fly histograms for switch vs background quality."""

    switch_gq: list[int] = field(default_factory=list)
    switch_dp: list[int] = field(default_factory=list)
    switch_rnc: list[str] = field(default_factory=list)
    bg_gq: list[int] = field(default_factory=list)
    bg_dp: list[int] = field(default_factory=list)
    bg_rnc: list[str] = field(default_factory=list)
    switch_low_gq: int = 0
    switch_low_dp: int = 0
    switch_rnc_i: int = 0
    switch_n: int = 0
    switch_sentinel: int = 0
    bg_low_gq: int = 0
    bg_low_dp: int = 0
    bg_rnc_i: int = 0
    bg_n: int = 0
    bg_sentinel: int = 0
    sites_scanned: int = 0
    sample_sites: int = 0
    switch_by_gt: dict[str, GtBin] = field(
        default_factory=lambda: {k: GtBin() for k in GT_CLASSES}
    )
    bg_by_gt: dict[str, GtBin] = field(
        default_factory=lambda: {k: GtBin() for k in GT_CLASSES}
    )


def parse_int(tok: str) -> Optional[int]:
    if tok in ("", ".", "./."):
        return None
    try:
        return int(tok)
    except ValueError:
        return None


def normalize_rnc(tok: str) -> str:
    if not tok or tok in (".", "./.", ".|."):
        return "."
    return tok


def rnc_has_I(rnc: str) -> bool:
    """True if GLnexus marked either allele incomplete (I)."""
    return "I" in (rnc or "")


def classify_gt(gt: str) -> str:
    """Map a VCF GT string to hom_ref / het / hom_alt / missing / other."""
    if not gt or gt in (".", "./.", ".|.", "./", "|"):
        return "missing"
    alleles = [a for a in gt.replace("|", "/").split("/") if a != ""]
    if not alleles or all(a == "." for a in alleles):
        return "missing"
    nums: list[int] = []
    for a in alleles:
        if a == ".":
            return "missing"
        try:
            nums.append(int(a))
        except ValueError:
            return "other"
    if all(n == 0 for n in nums):
        return "hom_ref"
    if len(set(nums)) == 1:
        return "hom_alt"
    return "het"


def is_glnexus_homref_sentinel(
    gt_class: str, gq: Optional[int], dp: Optional[int]
) -> bool:
    """GLnexus placeholder for cohort 0/0: DP=0 and GQ≤1."""
    if gt_class != "hom_ref":
        return False
    if dp not in (None, 0):
        return False
    if gq is None:
        return True
    return gq <= 1


def open_text(path: Path, mode: str = "rt"):
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def bcftools_bin() -> str:
    exe = shutil.which("bcftools")
    if not exe:
        raise RuntimeError("bcftools is required on PATH")
    return exe


def load_sample_ids(path: Optional[Path], vcf: Path) -> list[str]:
    if path is not None:
        ids = []
        with open(path) as fh:
            for line in fh:
                s = line.strip()
                if s and not s.startswith("#"):
                    ids.append(s.split()[0])
        if not ids:
            raise SystemExit(f"no sample IDs in {path}")
        return ids
    out = subprocess.check_output([bcftools_bin(), "query", "-l", str(vcf)], text=True)
    ids = [ln for ln in out.splitlines() if ln]
    if not ids:
        raise SystemExit(f"no samples in {vcf}")
    return ids


FORMAT_ID_RE = re.compile(r"^##FORMAT=<ID=([^,>]+)")


def vcf_format_ids(vcf: Path) -> set[str]:
    """FORMAT IDs declared on ##FORMAT= header lines (not substring hits)."""
    hdr = subprocess.check_output(
        [bcftools_bin(), "view", "-h", str(vcf)],
        text=True,
    )
    ids: set[str] = set()
    for line in hdr.splitlines():
        m = FORMAT_ID_RE.match(line)
        if m:
            ids.add(m.group(1))
    return ids


def vcf_has_format_tag(vcf: Path, tag: str) -> bool:
    return tag in vcf_format_ids(vcf)


# bcftools query errors if -f names a FORMAT tag absent from the header.
SAMPLE_FORMAT_ORDER = ("GQ", "DP", "AN1", "AN2", "RNC", "GT")


def sample_format_tags(present: dict[str, bool], *, an_only: bool = False) -> list[str]:
    """FORMAT tags to request, in query order. AN1/AN2 are required."""
    missing = [t for t in ("AN1", "AN2") if not present.get(t)]
    if missing:
        raise RuntimeError("VCF is missing FORMAT/" + " and ".join(missing))
    if an_only:
        return ["AN1", "AN2"]
    return [t for t in SAMPLE_FORMAT_ORDER if present.get(t)]


def query_format_string(tags: list[str], *, an_only: bool = False) -> str:
    """bcftools -f string. AN-only omits REF/ALT (unused for rate metrics)."""
    if an_only or tags == ["AN1", "AN2"]:
        return r"%CHROM\t%POS[\t%AN1\t%AN2]\n"
    inner = r"\t".join(f"%{t}" for t in tags)
    return r"%CHROM\t%POS\t%REF\t%ALT[\t" + inner + r"]\n"


def parse_sample_fields(
    parts: list[str],
    base: int,
    tags: list[str],
) -> tuple[Optional[int], Optional[int], Optional[int], Optional[int], str, str]:
    """Parse one sample's query fields into (gq, dp, an1, an2, rnc, gt)."""
    by_tag = {tag: parts[base + i] for i, tag in enumerate(tags)}
    gq = parse_int(by_tag["GQ"]) if "GQ" in by_tag else None
    dp = parse_int(by_tag["DP"]) if "DP" in by_tag else None
    an1 = parse_int(by_tag["AN1"])
    an2 = parse_int(by_tag["AN2"])
    rnc = normalize_rnc(by_tag["RNC"]) if "RNC" in by_tag else "."
    gt = by_tag["GT"] if "GT" in by_tag and by_tag["GT"] else "."
    return gq, dp, an1, an2, rnc, gt


def split_sample_chunks(samples: list[str], jobs: int) -> list[list[str]]:
    """Split samples into up to ``jobs`` contiguous non-empty chunks."""
    jobs = max(1, min(int(jobs), len(samples) or 1))
    if jobs == 1 or not samples:
        return [samples]
    base, rem = divmod(len(samples), jobs)
    chunks: list[list[str]] = []
    i = 0
    for c in range(jobs):
        sz = base + (1 if c < rem else 0)
        if sz:
            chunks.append(samples[i : i + sz])
            i += sz
    return chunks


def het_admixture_factor(props: Optional[Sequence[float]]) -> Optional[float]:
    """Return ``1 - sum p_k^2`` (probability two random draws differ in ancestry)."""
    if props is None:
        return None
    vals = [float(x) for x in props]
    if not vals:
        return None
    s = sum(vals)
    if s <= 0:
        return None
    # Tolerate un-normalized vectors (normalize).
    vals = [x / s for x in vals]
    return 1.0 - sum(x * x for x in vals)


def expected_mean_tract_mb(t_gen: float, rec_cm_per_mb: float = DEFAULT_REC_CM_PER_MB) -> Optional[float]:
    """Mean admixture-tract length in Mb for generations-since-admixture T.

    This is the unconditional tract length scale ``100 / (T * rec)``. The
    expected *switch* rate also multiplies by ``(1 - sum p_k^2)``.
    """
    if t_gen is None or t_gen <= 0 or rec_cm_per_mb <= 0:
        return None
    return 100.0 / (t_gen * rec_cm_per_mb)


def expected_switches_per_hap_per_mb(
    t_gen: float,
    props: Optional[Sequence[float]] = None,
    rec_cm_per_mb: float = DEFAULT_REC_CM_PER_MB,
) -> Optional[float]:
    """Expected ancestry switches per haplotype per Mb.

    Rate = ``(T * rec / 100) * (1 - sum p_k^2)``. When ``props`` is omitted,
    the heterozygosity factor is treated as 1.0 (upper bound / legacy).
    ``rate_over_expected`` is monotone in T; it falsifies, it does not select T.
    """
    if t_gen is None or t_gen <= 0 or rec_cm_per_mb <= 0:
        return None
    base = (t_gen * rec_cm_per_mb) / 100.0
    if props is None:
        return base
    factor = het_admixture_factor(props)
    if factor is None:
        return None
    return base * factor


def implied_t_gen(
    switches_per_hap_per_mb: Optional[float],
    rec_cm_per_mb: float = DEFAULT_REC_CM_PER_MB,
) -> Optional[float]:
    """Legacy implied T assuming heterozygosity factor 1 (no props). Prefer
    :func:`implied_T_given_props` when proportions are known."""
    if switches_per_hap_per_mb is None or switches_per_hap_per_mb <= 0 or rec_cm_per_mb <= 0:
        return None
    return switches_per_hap_per_mb * 100.0 / rec_cm_per_mb


def implied_T_given_props(
    switches_per_hap_per_mb: Optional[float],
    props: Optional[Sequence[float]],
    rec_cm_per_mb: float = DEFAULT_REC_CM_PER_MB,
) -> Optional[float]:
    """Invert ``rate = (T * rec / 100) * (1 - sum p^2)`` for T."""
    if switches_per_hap_per_mb is None or switches_per_hap_per_mb <= 0 or rec_cm_per_mb <= 0:
        return None
    factor = het_admixture_factor(props)
    if factor is None or factor <= 0:
        return None
    return switches_per_hap_per_mb * 100.0 / (rec_cm_per_mb * factor)


def span_ok_for_t(
    span_mb: Optional[float],
    t_gen: Optional[float],
    *,
    min_mean_tracts: float = 3.0,
    rec_cm_per_mb: float = DEFAULT_REC_CM_PER_MB,
) -> Optional[bool]:
    """True when ``span_mb >= min_mean_tracts * (100 / T)`` (mean tract lengths)."""
    if span_mb is None or t_gen is None or t_gen <= 0:
        return None
    mean_mb = expected_mean_tract_mb(t_gen, rec_cm_per_mb)
    if mean_mb is None:
        return None
    return span_mb >= min_mean_tracts * mean_mb


def iter_annotated_sites(
    vcf: Path,
    samples: list[str],
    region: Optional[str] = None,
    *,
    an_only: bool = False,
) -> Iterator[tuple[str, int, str, str, list[tuple[Optional[int], Optional[int], Optional[int], Optional[int], str, str]]]]:
    """Yield (chrom, pos, ref, alt, per_sample[(gq, dp, an1, an2, rnc, gt), ...]).

    With ``an_only=True``, REF/ALT are \".\" and GQ/DP/RNC/GT are unused placeholders.
    """
    format_ids = vcf_format_ids(vcf)
    present = {tag: tag in format_ids for tag in SAMPLE_FORMAT_ORDER}
    tags = sample_format_tags(present, an_only=an_only)
    if an_only:
        print(
            "GQ/DP unavailable: --an-only mode (FORMAT GQ/DP not queried). "
            "gq_dp_available=false; quality histograms will be null, not zeros.",
            file=sys.stderr,
            flush=True,
        )
    elif not present.get("GQ") or not present.get("DP"):
        print(
            "GQ/DP unavailable: FORMAT/GQ or DP missing from VCF. "
            "gq_dp_available=false; quality histograms will be null, not zeros. "
            "Switch rates and tracts still run on AN1/AN2.",
            file=sys.stderr,
            flush=True,
        )
    if not an_only and not present.get("RNC"):
        print(
            "warning: FORMAT/RNC missing from VCF header; RNC-I stats will be empty. "
            "Use AnnotateFlareGqDp with format_tags=GQ,DP,RNC and a fresh local copy "
            "(the notebook caches RUN_DIR/<vcf-name>).",
            file=sys.stderr,
            flush=True,
        )
    fmt = query_format_string(tags, an_only=an_only)
    fields_per = 2 if an_only else len(tags)
    print(f"bcftools query FORMAT: {','.join(tags)}" + (" (an-only)" if an_only else ""), flush=True)
    cmd = [
        bcftools_bin(),
        "query",
        "-f",
        fmt,
        "-s",
        ",".join(samples),
    ]
    if region:
        cmd.extend(["-r", region])
    cmd.append(str(vcf))
    proc = subprocess.Popen(
        cmd,
        stdout=subprocess.PIPE,
        stderr=sys.stderr,
        text=True,
        bufsize=1 << 20,
    )
    assert proc.stdout is not None
    n_samp = len(samples)
    try:
        for line_no, line in enumerate(proc.stdout, 1):
            parts = line.rstrip("\n").split("\t")
            if an_only:
                expect = 2 + fields_per * n_samp
                if len(parts) < expect:
                    raise RuntimeError(
                        f"line {line_no}: expected {expect} fields, got {len(parts)}. "
                        "Need FORMAT AN1/AN2."
                    )
                chrom, pos_s = parts[0], parts[1]
                pos = int(pos_s)
                per = []
                base = 2
                for _ in range(n_samp):
                    an1 = parse_int(parts[base])
                    an2 = parse_int(parts[base + 1])
                    per.append((None, None, an1, an2, ".", "."))
                    base += 2
                yield chrom, pos, ".", ".", per
            else:
                expect = 4 + fields_per * n_samp
                if len(parts) < expect:
                    raise RuntimeError(
                        f"line {line_no}: expected {expect} fields, got {len(parts)}. "
                        "Need FORMAT AN1/AN2 (GQ/DP/RNC/GT optional; omitted from -f if absent)."
                    )
                chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
                pos = int(pos_s)
                per = []
                base = 4
                for _ in range(n_samp):
                    per.append(parse_sample_fields(parts, base, tags))
                    base += fields_per
                yield chrom, pos, ref, alt, per
    finally:
        proc.stdout.close()
        rc = proc.wait()
        if rc not in (0, None):
            raise RuntimeError(f"bcftools query failed (exit {rc})")


def mark_flickers(switches: list[SwitchEvent], max_span_bp: int) -> None:
    """Flag A→B→A (or B→A→B) on the same hap within max_span_bp."""
    by_hap: dict[int, list[SwitchEvent]] = {1: [], 2: []}
    for ev in switches:
        by_hap[ev.hap].append(ev)
    for hap_events in by_hap.values():
        for i in range(len(hap_events) - 1):
            a = hap_events[i]
            b = hap_events[i + 1]
            if a.chrom != b.chrom:
                continue
            span = b.pos - a.pos
            if span <= 0 or span > max_span_bp:
                continue
            # Reversion: to of first == from of second and from of first == to of second
            if a.anc_to == b.anc_from and a.anc_from == b.anc_to:
                a.is_flicker = True
                b.is_flicker = True
                a.flicker_span_bp = span
                b.flicker_span_bp = span


def maybe_record_switch(
    tracker: SampleTracker,
    sample: str,
    hap: int,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    anc: Optional[int],
    gq: Optional[int],
    dp: Optional[int],
    rnc: str,
    gt: str = ".",
    gt_class: str = "missing",
    is_sentinel: bool = False,
) -> Optional[SwitchEvent]:
    state = tracker.hap1 if hap == 1 else tracker.hap2
    if anc is None:
        return None
    ev = None
    if state.anc is not None and state.anc != anc and state.pos is not None:
        prev_het_gap = (
            None if tracker.last_het_pos is None else pos - tracker.last_het_pos
        )
        ev = SwitchEvent(
            sample=sample,
            hap=hap,
            chrom=chrom,
            pos=pos,
            ref=ref,
            alt=alt,
            anc_from=state.anc,
            anc_to=anc,
            gq=gq,
            dp=dp,
            rnc=rnc,
            prev_pos=state.pos,
            bp_gap=pos - state.pos,
            gt=gt,
            gt_class=gt_class,
            is_glnexus_homref_sentinel=is_sentinel,
            prev_het_pos=tracker.last_het_pos,
            prev_het_gq=tracker.last_het_gq,
            prev_het_dp=tracker.last_het_dp,
            prev_het_gap=prev_het_gap,
            prev_het_rnc=tracker.last_het_rnc,
        )
        tracker.switches.append(ev)
        tracker.pending_next_het.append(ev)
    state.anc = anc
    state.pos = pos
    state.gq = gq
    state.dp = dp
    return ev


def maybe_record_switch_an(
    tracker: SampleTracker,
    sample: str,
    hap: int,
    chrom: str,
    pos: int,
    anc: Optional[int],
) -> Optional[SwitchEvent]:
    """AN-only switch recorder (no GQ/DP/het bookkeeping)."""
    state = tracker.hap1 if hap == 1 else tracker.hap2
    if anc is None:
        return None
    ev = None
    if state.anc is not None and state.anc != anc and state.pos is not None:
        ev = SwitchEvent(
            sample=sample,
            hap=hap,
            chrom=chrom,
            pos=pos,
            ref=".",
            alt=".",
            anc_from=state.anc,
            anc_to=anc,
            gq=None,
            dp=None,
            rnc=".",
            prev_pos=state.pos,
            bp_gap=pos - state.pos,
        )
        tracker.switches.append(ev)
    state.anc = anc
    state.pos = pos
    return ev


def reset_tracker_chrom(tracker: SampleTracker) -> None:
    """Drop hap ancestry and flanking-het state across chromosome boundaries."""
    tracker.hap1 = HapState()
    tracker.hap2 = HapState()
    tracker.last_het_pos = None
    tracker.last_het_gq = None
    tracker.last_het_dp = None
    tracker.last_het_rnc = "."
    tracker.pending_next_het = []


def note_het(
    tracker: SampleTracker,
    chrom: str,
    pos: int,
    gq: Optional[int],
    dp: Optional[int],
    rnc: str = ".",
) -> None:
    """Fill next-het on pending switches strictly before this site; then remember it."""
    still_pending: list[SwitchEvent] = []
    for ev in tracker.pending_next_het:
        if ev.chrom == chrom and ev.pos < pos and ev.next_het_pos is None:
            ev.next_het_pos = pos
            ev.next_het_gq = gq
            ev.next_het_dp = dp
            ev.next_het_gap = pos - ev.pos
            ev.next_het_rnc = rnc
        else:
            still_pending.append(ev)
    tracker.pending_next_het = still_pending
    tracker.last_het_pos = pos
    tracker.last_het_gq = gq
    tracker.last_het_dp = dp
    tracker.last_het_rnc = rnc


def assign_nearest_het(ev: SwitchEvent) -> None:
    """Nearest informative het: the switch site if it is a het, else closer flanking het."""
    if ev.gt_class == "het":
        ev.nearest_het_pos = ev.pos
        ev.nearest_het_gq = ev.gq
        ev.nearest_het_dp = ev.dp
        ev.nearest_het_gap = 0
        ev.nearest_het_side = "current"
        ev.nearest_het_rnc = ev.rnc
        return
    candidates: list[tuple[str, int, Optional[int], Optional[int], int, str]] = []
    if ev.prev_het_pos is not None:
        candidates.append(
            (
                "prev",
                ev.prev_het_pos,
                ev.prev_het_gq,
                ev.prev_het_dp,
                ev.pos - ev.prev_het_pos,
                ev.prev_het_rnc,
            )
        )
    if ev.next_het_pos is not None:
        candidates.append(
            (
                "next",
                ev.next_het_pos,
                ev.next_het_gq,
                ev.next_het_dp,
                ev.next_het_pos - ev.pos,
                ev.next_het_rnc,
            )
        )
    if not candidates:
        ev.nearest_het_side = "none"
        return
    side, het_pos, het_gq, het_dp, gap, het_rnc = min(
        candidates, key=lambda c: (c[4], 0 if c[0] == "prev" else 1)
    )
    ev.nearest_het_side = side
    ev.nearest_het_pos = het_pos
    ev.nearest_het_gq = het_gq
    ev.nearest_het_dp = het_dp
    ev.nearest_het_gap = gap
    ev.nearest_het_rnc = het_rnc


def process_sample_site(
    tracker: SampleTracker,
    sample: str,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    an1: Optional[int],
    an2: Optional[int],
    gq: Optional[int],
    dp: Optional[int],
    rnc: str,
    gt: str,
    *,
    an_only: bool = False,
) -> tuple[Optional[SwitchEvent], Optional[SwitchEvent], str, bool]:
    """Record switches, then update last/next het. Returns (ev1, ev2, gt_class, sentinel)."""
    if tracker.last_chrom is not None and tracker.last_chrom != chrom:
        reset_tracker_chrom(tracker)
    tracker.last_chrom = chrom
    if an_only:
        ev1 = maybe_record_switch_an(tracker, sample, 1, chrom, pos, an1)
        ev2 = maybe_record_switch_an(tracker, sample, 2, chrom, pos, an2)
        return ev1, ev2, "missing", False
    gt_class = classify_gt(gt)
    sentinel = is_glnexus_homref_sentinel(gt_class, gq, dp)
    ev1 = maybe_record_switch(
        tracker, sample, 1, chrom, pos, ref, alt, an1, gq, dp, rnc,
        gt=gt, gt_class=gt_class, is_sentinel=sentinel,
    )
    ev2 = maybe_record_switch(
        tracker, sample, 2, chrom, pos, ref, alt, an2, gq, dp, rnc,
        gt=gt, gt_class=gt_class, is_sentinel=sentinel,
    )
    if gt_class == "het":
        note_het(tracker, chrom, pos, gq, dp, rnc)
    return ev1, ev2, gt_class, sentinel


def update_qc(
    qc: QcAccum,
    gq: Optional[int],
    dp: Optional[int],
    rnc: str,
    *,
    is_switch: bool,
    gq_thr: int,
    dp_thr: int,
    bg_keep_every: int,
    bg_counter: list[int],
    gt_class: str = "missing",
    is_sentinel: bool = False,
) -> None:
    def low_gq(v: Optional[int]) -> bool:
        return v is not None and v < gq_thr

    def low_dp(v: Optional[int]) -> bool:
        return v is not None and v < dp_thr

    def fill_bin(bin_: GtBin, keep_hist: bool) -> None:
        bin_.n += 1
        if low_gq(gq):
            bin_.low_gq += 1
        if low_dp(dp):
            bin_.low_dp += 1
        if is_sentinel:
            bin_.sentinel += 1
        if rnc_has_I(rnc):
            bin_.rnc_i += 1
        if keep_hist:
            if gq is not None and len(bin_.gq) < 2_000_000:
                bin_.gq.append(gq)
            if dp is not None and len(bin_.dp) < 2_000_000:
                bin_.dp.append(dp)

    rnc_bad = rnc_has_I(rnc)
    klass = gt_class if gt_class in qc.switch_by_gt else "other"
    if is_switch:
        qc.switch_n += 1
        if is_sentinel:
            qc.switch_sentinel += 1
        if gq is not None:
            qc.switch_gq.append(gq)
        if dp is not None:
            qc.switch_dp.append(dp)
        if len(qc.switch_rnc) < 2_000_000:
            qc.switch_rnc.append(normalize_rnc(rnc))
        if low_gq(gq):
            qc.switch_low_gq += 1
        if low_dp(dp):
            qc.switch_low_dp += 1
        if rnc_bad:
            qc.switch_rnc_i += 1
        fill_bin(qc.switch_by_gt[klass], keep_hist=True)
        return

    bg_counter[0] += 1
    if bg_counter[0] % bg_keep_every != 0:
        return
    qc.bg_n += 1
    if is_sentinel:
        qc.bg_sentinel += 1
    if gq is not None:
        if len(qc.bg_gq) < 2_000_000:
            qc.bg_gq.append(gq)
    if dp is not None:
        if len(qc.bg_dp) < 2_000_000:
            qc.bg_dp.append(dp)
    if len(qc.bg_rnc) < 2_000_000:
        qc.bg_rnc.append(normalize_rnc(rnc))
    if low_gq(gq):
        qc.bg_low_gq += 1
    if low_dp(dp):
        qc.bg_low_dp += 1
    if rnc_bad:
        qc.bg_rnc_i += 1
    fill_bin(qc.bg_by_gt[klass], keep_hist=True)


def mean(xs: list[int]) -> Optional[float]:
    if not xs:
        return None
    return sum(xs) / len(xs)


def percentile(xs: list[int], p: float) -> Optional[float]:
    if not xs:
        return None
    ys = sorted(xs)
    if len(ys) == 1:
        return float(ys[0])
    k = (len(ys) - 1) * p
    f = math.floor(k)
    c = math.ceil(k)
    if f == c:
        return float(ys[int(k)])
    return ys[f] * (c - k) + ys[c] * (k - f)


def build_tracts(
    events: list[SwitchEvent],
    samples: list[str],
    last_anc: dict[tuple[str, int], Optional[int]],
    pos_min: dict[str, int],
    pos_max: dict[str, int],
) -> list[AncTract]:
    """Window-censored ancestry stretches, including haps with no switch."""
    by: dict[tuple[str, int, str], list[SwitchEvent]] = defaultdict(list)
    for ev in events:
        by[(ev.sample, ev.hap, ev.chrom)].append(ev)
    chroms = [c for c in pos_min if c in pos_max and pos_max[c] > pos_min[c]]
    tracts: list[AncTract] = []
    for sample in samples:
        for hap in (1, 2):
            for chrom in chroms:
                lo = pos_min[chrom]
                hi = pos_max[chrom]
                evs = sorted(by.get((sample, hap, chrom), []), key=lambda e: e.pos)
                if not evs:
                    anc = last_anc.get((sample, hap))
                    if anc is None:
                        continue
                    tracts.append(
                        AncTract(
                            sample=sample,
                            hap=hap,
                            chrom=chrom,
                            start=lo,
                            end=hi,
                            anc=anc,
                            left_censored=True,
                            right_censored=True,
                            n_switches_on_hap=0,
                        )
                    )
                    continue
                n_sw = len(evs)
                first_end = evs[0].pos
                if first_end > lo:
                    tracts.append(
                        AncTract(
                            sample=sample,
                            hap=hap,
                            chrom=chrom,
                            start=lo,
                            end=first_end,
                            anc=evs[0].anc_from,
                            left_censored=True,
                            right_censored=False,
                            n_switches_on_hap=n_sw,
                        )
                    )
                for i in range(len(evs) - 1):
                    a = evs[i]
                    b = evs[i + 1]
                    if b.pos <= a.pos:
                        continue
                    tracts.append(
                        AncTract(
                            sample=sample,
                            hap=hap,
                            chrom=chrom,
                            start=a.pos,
                            end=b.pos,
                            anc=a.anc_to,
                            left_censored=False,
                            right_censored=False,
                            n_switches_on_hap=n_sw,
                            is_flicker_interval=bool(a.is_flicker and b.is_flicker),
                        )
                    )
                last = evs[-1]
                if hi > last.pos:
                    tracts.append(
                        AncTract(
                            sample=sample,
                            hap=hap,
                            chrom=chrom,
                            start=last.pos,
                            end=hi,
                            anc=last.anc_to,
                            left_censored=False,
                            right_censored=True,
                            n_switches_on_hap=n_sw,
                        )
                    )
    return tracts


def scan_vcf(
    vcf: Path,
    samples: list[str],
    *,
    gq_threshold: int,
    dp_threshold: int,
    flicker_max_bp: int,
    bg_keep_every: int,
    max_sites: Optional[int] = None,
    region: Optional[str] = None,
    an_only: bool = False,
    rates_only: bool = False,
) -> tuple[list[SwitchEvent], QcAccum, dict, list[AncTract]]:
    trackers = {s: SampleTracker() for s in samples}
    qc = QcAccum()
    bg_counter = [0]
    chrom_sites: dict[str, int] = defaultdict(int)
    pos_min: dict[str, int] = {}
    pos_max: dict[str, int] = {}
    skip_qc = an_only or rates_only

    for chrom, pos, ref, alt, per in iter_annotated_sites(
        vcf, samples, region=region, an_only=an_only
    ):
        qc.sites_scanned += 1
        chrom_sites[chrom] += 1
        if chrom not in pos_min:
            pos_min[chrom] = pos
        pos_max[chrom] = pos
        for sample, (gq, dp, an1, an2, rnc, gt) in zip(samples, per):
            tr = trackers[sample]
            qc.sample_sites += 1
            ev1, ev2, gt_class, sentinel = process_sample_site(
                tr, sample, chrom, pos, ref, alt, an1, an2, gq, dp, rnc, gt,
                an_only=an_only,
            )
            if skip_qc:
                continue
            is_switch = ev1 is not None or ev2 is not None
            # Quality is genotype-level; attribute once per sample-site even if both haps switch.
            update_qc(
                qc,
                gq,
                dp,
                rnc,
                is_switch=is_switch,
                gq_thr=gq_threshold,
                dp_thr=dp_threshold,
                bg_keep_every=bg_keep_every,
                bg_counter=bg_counter,
                gt_class=gt_class,
                is_sentinel=sentinel,
            )
        if max_sites is not None and qc.sites_scanned >= max_sites:
            break

    all_switches: list[SwitchEvent] = []
    for tr in trackers.values():
        if not an_only:
            for ev in tr.switches:
                assign_nearest_het(ev)
        mark_flickers(tr.switches, flicker_max_bp)
        all_switches.extend(tr.switches)

    last_anc: dict[tuple[str, int], Optional[int]] = {}
    for sample, tr in trackers.items():
        last_anc[(sample, 1)] = tr.hap1.anc
        last_anc[(sample, 2)] = tr.hap2.anc
    tracts = build_tracts(all_switches, samples, last_anc, pos_min, pos_max)
    span_bp = sum(pos_max[c] - pos_min[c] for c in pos_min if c in pos_max)

    format_ids = vcf_format_ids(vcf)
    gq_dp_available = (not an_only) and ("GQ" in format_ids) and ("DP" in format_ids)

    meta = {
        "n_samples": len(samples),
        "sites_scanned": qc.sites_scanned,
        "chrom_sites": dict(chrom_sites),
        "n_switches": len(all_switches),
        "n_flicker_switches": sum(1 for e in all_switches if e.is_flicker),
        "gq_threshold": gq_threshold,
        "dp_threshold": dp_threshold,
        "flicker_max_bp": flicker_max_bp,
        "bg_keep_every": bg_keep_every,
        "region": region,
        "pos_min": dict(pos_min),
        "pos_max": dict(pos_max),
        "span_bp": span_bp,
        "an_only": an_only,
        "rates_only": rates_only,
        "gq_dp_available": gq_dp_available,
    }
    return all_switches, qc, meta, tracts


def merge_qc(parts: list[QcAccum]) -> QcAccum:
    """Concatenate per-chunk QC accumulators (sites_scanned taken as max)."""
    out = QcAccum()
    if not parts:
        return out
    out.sites_scanned = max(q.sites_scanned for q in parts)
    for q in parts:
        out.sample_sites += q.sample_sites
        out.switch_n += q.switch_n
        out.switch_low_gq += q.switch_low_gq
        out.switch_low_dp += q.switch_low_dp
        out.switch_rnc_i += q.switch_rnc_i
        out.switch_sentinel += q.switch_sentinel
        out.bg_n += q.bg_n
        out.bg_low_gq += q.bg_low_gq
        out.bg_low_dp += q.bg_low_dp
        out.bg_rnc_i += q.bg_rnc_i
        out.bg_sentinel += q.bg_sentinel
        out.switch_gq.extend(q.switch_gq)
        out.switch_dp.extend(q.switch_dp)
        out.switch_rnc.extend(q.switch_rnc)
        out.bg_gq.extend(q.bg_gq)
        out.bg_dp.extend(q.bg_dp)
        out.bg_rnc.extend(q.bg_rnc)
        for klass in GT_CLASSES:
            sw = out.switch_by_gt.setdefault(klass, GtBin())
            bg = out.bg_by_gt.setdefault(klass, GtBin())
            sw_q = q.switch_by_gt.get(klass, GtBin())
            bg_q = q.bg_by_gt.get(klass, GtBin())
            for dst, src in ((sw, sw_q), (bg, bg_q)):
                dst.n += src.n
                dst.low_gq += src.low_gq
                dst.low_dp += src.low_dp
                dst.sentinel += src.sentinel
                dst.rnc_i += src.rnc_i
                dst.gq.extend(src.gq)
                dst.dp.extend(src.dp)
    return out


def merge_scan_results(
    parts: list[tuple[list[SwitchEvent], QcAccum, dict, list[AncTract]]],
) -> tuple[list[SwitchEvent], QcAccum, dict, list[AncTract]]:
    if len(parts) == 1:
        return parts[0]
    events: list[SwitchEvent] = []
    tracts: list[AncTract] = []
    qcs: list[QcAccum] = []
    metas: list[dict] = []
    for ev, qc, meta, tr in parts:
        events.extend(ev)
        tracts.extend(tr)
        qcs.append(qc)
        metas.append(meta)
    meta0 = dict(metas[0])
    meta0["n_samples"] = sum(m["n_samples"] for m in metas)
    meta0["n_switches"] = sum(m["n_switches"] for m in metas)
    meta0["n_flicker_switches"] = sum(m["n_flicker_switches"] for m in metas)
    meta0["sites_scanned"] = max(m["sites_scanned"] for m in metas)
    meta0["jobs"] = len(parts)
    qc = merge_qc(qcs)
    # rates-only workers leave QcAccum.sites_scanned=0; trust meta.
    qc.sites_scanned = meta0["sites_scanned"]
    return events, qc, meta0, tracts


def _scan_vcf_worker(payload: dict) -> tuple[list[SwitchEvent], QcAccum, dict, list[AncTract]]:
    """Picklable worker for ProcessPoolExecutor."""
    return scan_vcf(
        Path(payload["vcf"]),
        payload["samples"],
        gq_threshold=payload["gq_threshold"],
        dp_threshold=payload["dp_threshold"],
        flicker_max_bp=payload["flicker_max_bp"],
        bg_keep_every=payload["bg_keep_every"],
        max_sites=payload.get("max_sites"),
        region=payload.get("region"),
        an_only=payload.get("an_only", False),
        rates_only=payload.get("rates_only", False),
    )


def scan_vcf_parallel(
    vcf: Path,
    samples: list[str],
    *,
    jobs: int,
    gq_threshold: int,
    dp_threshold: int,
    flicker_max_bp: int,
    bg_keep_every: int,
    max_sites: Optional[int] = None,
    region: Optional[str] = None,
    an_only: bool = False,
    rates_only: bool = False,
) -> tuple[list[SwitchEvent], QcAccum, dict, list[AncTract]]:
    chunks = split_sample_chunks(samples, jobs)
    if len(chunks) == 1:
        return scan_vcf(
            vcf,
            samples,
            gq_threshold=gq_threshold,
            dp_threshold=dp_threshold,
            flicker_max_bp=flicker_max_bp,
            bg_keep_every=bg_keep_every,
            max_sites=max_sites,
            region=region,
            an_only=an_only,
            rates_only=rates_only,
        )
    print(
        f"Parallel scan: {len(chunks)} jobs × ~{len(chunks[0])} samples "
        f"(of {len(samples)})",
        flush=True,
    )
    payloads = [
        {
            "vcf": str(vcf),
            "samples": chunk,
            "gq_threshold": gq_threshold,
            "dp_threshold": dp_threshold,
            "flicker_max_bp": flicker_max_bp,
            "bg_keep_every": bg_keep_every,
            "max_sites": max_sites,
            "region": region,
            "an_only": an_only,
            "rates_only": rates_only,
        }
        for chunk in chunks
    ]
    with ProcessPoolExecutor(max_workers=len(chunks)) as pool:
        parts = list(pool.map(_scan_vcf_worker, payloads))
    return merge_scan_results(parts)


def write_switches_tsv(path: Path, events: list[SwitchEvent]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        def ni(v: Optional[int]) -> str:
            return "" if v is None else str(v)

        fh.write(
            "sample\thap\tchrom\tpos\tref\talt\tgt\tgt_class\t"
            "is_glnexus_homref_sentinel\tanc_from\tanc_to\t"
            "anc_from_label\tanc_to_label\tgq\tdp\trnc\tprev_pos\tbp_gap\t"
            "is_flicker\tflicker_span_bp\t"
            "prev_het_pos\tprev_het_gq\tprev_het_dp\tprev_het_gap\tprev_het_rnc\t"
            "next_het_pos\tnext_het_gq\tnext_het_dp\tnext_het_gap\tnext_het_rnc\t"
            "nearest_het_pos\tnearest_het_gq\tnearest_het_dp\tnearest_het_gap\t"
            "nearest_het_side\tnearest_het_rnc\n"
        )
        for e in events:
            fh.write(
                f"{e.sample}\t{e.hap}\t{e.chrom}\t{e.pos}\t{e.ref}\t{e.alt}\t"
                f"{e.gt}\t{e.gt_class}\t{int(e.is_glnexus_homref_sentinel)}\t"
                f"{e.anc_from}\t{e.anc_to}\t"
                f"{ANCESTRY.get(e.anc_from, '?')}\t{ANCESTRY.get(e.anc_to, '?')}\t"
                f"{ni(e.gq)}\t{ni(e.dp)}\t"
                f"{e.rnc}\t{e.prev_pos}\t{e.bp_gap}\t"
                f"{int(e.is_flicker)}\t"
                f"{ni(e.flicker_span_bp)}\t"
                f"{ni(e.prev_het_pos)}\t{ni(e.prev_het_gq)}\t{ni(e.prev_het_dp)}\t"
                f"{ni(e.prev_het_gap)}\t{e.prev_het_rnc}\t"
                f"{ni(e.next_het_pos)}\t{ni(e.next_het_gq)}\t{ni(e.next_het_dp)}\t"
                f"{ni(e.next_het_gap)}\t{e.next_het_rnc}\t"
                f"{ni(e.nearest_het_pos)}\t{ni(e.nearest_het_gq)}\t"
                f"{ni(e.nearest_het_dp)}\t{ni(e.nearest_het_gap)}\t"
                f"{e.nearest_het_side}\t{e.nearest_het_rnc}\n"
            )


def write_switches_slim_tsv(path: Path, events: list[SwitchEvent]) -> None:
    """Compact switch table for by-pop rate rollups (sample + flicker only)."""
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        fh.write("sample\thap\tis_flicker\n")
        for e in events:
            fh.write(f"{e.sample}\t{e.hap}\t{int(e.is_flicker)}\n")


def write_tracts_tsv(path: Path, tracts: list[AncTract]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        fh.write(
            "sample\thap\tchrom\tstart\tend\tlength_bp\tanc\tanc_label\t"
            "left_censored\tright_censored\tn_switches_on_hap\tis_flicker_interval\n"
        )
        for t in tracts:
            anc_s = "" if t.anc is None else str(t.anc)
            label = ANCESTRY.get(t.anc, "?") if t.anc is not None else ""
            fh.write(
                f"{t.sample}\t{t.hap}\t{t.chrom}\t{t.start}\t{t.end}\t{t.length_bp}\t"
                f"{anc_s}\t{label}\t{int(t.left_censored)}\t{int(t.right_censored)}\t"
                f"{t.n_switches_on_hap}\t{int(t.is_flicker_interval)}\n"
            )


def tract_stats(
    tracts: list[AncTract],
    meta: dict,
    events: list[SwitchEvent],
    rec_cm_per_mb: float = DEFAULT_REC_CM_PER_MB,
) -> dict:
    def frac(num: int, den: int) -> Optional[float]:
        if den == 0:
            return None
        return num / den

    n_samples = meta["n_samples"]
    n_haps = 2 * n_samples
    span_bp = meta.get("span_bp") or 0
    span_mb = span_bp / 1e6 if span_bp else None
    n_switches = len(events)
    n_flicker = sum(1 for e in events if e.is_flicker)
    n_haps_with_switch = len({(e.sample, e.hap) for e in events})
    rate = None
    if n_haps and span_mb:
        rate = n_switches / n_haps / span_mb
    flicker_rate = None
    if n_haps and span_mb:
        flicker_rate = n_flicker / n_haps / span_mb
    lengths = [t.length_bp for t in tracts if t.length_bp > 0]
    complete = [
        t.length_bp
        for t in tracts
        if t.length_bp > 0 and not t.left_censored and not t.right_censored
    ]
    flicker_lens = [t.length_bp for t in tracts if t.is_flicker_interval and t.length_bp > 0]
    return {
        "n_tracts": len(tracts),
        "n_complete_tracts": len(complete),
        "n_haps": n_haps,
        "n_haps_with_switch": n_haps_with_switch,
        "frac_haps_with_switch": frac(n_haps_with_switch, n_haps),
        "span_bp": span_bp,
        "span_mb": span_mb,
        "switches_per_hap_per_mb": rate,
        "flicker_frac": frac(n_flicker, n_switches),
        "flicker_per_hap_per_mb": flicker_rate,
        "implied_t_gen": implied_t_gen(rate, rec_cm_per_mb),
        "rec_cm_per_mb": rec_cm_per_mb,
        "median_length_bp": percentile(lengths, 0.5),
        "median_complete_length_bp": percentile(complete, 0.5),
        "p10_length_bp": percentile(lengths, 0.1),
        "p90_length_bp": percentile(lengths, 0.9),
        "median_flicker_interval_bp": percentile(flicker_lens, 0.5),
        "n_flicker_intervals": len(flicker_lens),
    }


def switch_fails_call_qc(
    ev: SwitchEvent,
    *,
    gq_threshold: int,
    dp_threshold: int,
    drop_rnc_i: bool = True,
) -> bool:
    """True if this switch should be dropped under Part 7 call-QC rules."""
    if ev.gq is not None and ev.gq < gq_threshold:
        return True
    if ev.dp is not None and ev.dp < dp_threshold:
        return True
    if drop_rnc_i and rnc_has_I(ev.rnc):
        return True
    return False


def filter_switches_call_qc(
    events: list[SwitchEvent],
    *,
    gq_threshold: int,
    dp_threshold: int,
    drop_rnc_i: bool = True,
) -> tuple[list[SwitchEvent], list[SwitchEvent]]:
    """Split switches into (pass, fail) under GQ/DP/RNC call-QC."""
    passed: list[SwitchEvent] = []
    failed: list[SwitchEvent] = []
    for ev in events:
        if switch_fails_call_qc(
            ev,
            gq_threshold=gq_threshold,
            dp_threshold=dp_threshold,
            drop_rnc_i=drop_rnc_i,
        ):
            failed.append(ev)
        else:
            passed.append(ev)
    return passed, failed


def call_qc_filtered_tract_metrics(
    events: list[SwitchEvent],
    meta: dict,
    *,
    gq_threshold: Optional[int] = None,
    dp_threshold: Optional[int] = None,
    drop_rnc_i: bool = True,
    rec_cm_per_mb: float = 1.0,
) -> dict:
    """Switch rates after dropping poor GQ/DP/RNC switch sites (Part 7.1).

    Does not rebuild tracts (ancestry path still includes those sites in FLARE
    output); this only re-scores annotation noise attributable to bad calls.
    """
    thr_gq = int(gq_threshold if gq_threshold is not None else meta.get("gq_threshold", 20))
    thr_dp = int(dp_threshold if dp_threshold is not None else meta.get("dp_threshold", 10))
    passed, failed = filter_switches_call_qc(
        events,
        gq_threshold=thr_gq,
        dp_threshold=thr_dp,
        drop_rnc_i=drop_rnc_i,
    )
    base = tract_stats([], meta, passed, rec_cm_per_mb=rec_cm_per_mb)
    n_all = len(events)
    n_fail = len(failed)
    n_pass = len(passed)

    def frac(num: int, den: int) -> Optional[float]:
        if den == 0:
            return None
        return num / den

    return {
        "gq_threshold": thr_gq,
        "dp_threshold": thr_dp,
        "drop_rnc_i": drop_rnc_i,
        "n_switches_all": n_all,
        "n_switches_fail_call_qc": n_fail,
        "n_switches_pass_call_qc": n_pass,
        "frac_switches_fail_call_qc": frac(n_fail, n_all),
        "switches_per_hap_per_mb": base.get("switches_per_hap_per_mb"),
        "flicker_frac": base.get("flicker_frac"),
        "flicker_per_hap_per_mb": base.get("flicker_per_hap_per_mb"),
        "implied_t_gen": base.get("implied_t_gen"),
        "n_haps_with_switch": base.get("n_haps_with_switch"),
        "frac_haps_with_switch": base.get("frac_haps_with_switch"),
    }


def summarize(
    qc: QcAccum,
    meta: dict,
    events: list[SwitchEvent],
    tracts: Optional[list[AncTract]] = None,
) -> dict:
    def frac(num: int, den: int) -> Optional[float]:
        if den == 0:
            return None
        return num / den

    flicker = [e for e in events if e.is_flicker]
    sustained = [e for e in events if not e.is_flicker]

    def qual_block(label: str, subset: list[SwitchEvent]) -> dict:
        gqs = [e.gq for e in subset if e.gq is not None]
        dps = [e.dp for e in subset if e.dp is not None]
        thr_gq = meta["gq_threshold"]
        thr_dp = meta["dp_threshold"]
        return {
            "n": len(subset),
            "mean_gq": mean(gqs),
            "median_gq": percentile(gqs, 0.5),
            "frac_gq_lt_thr": frac(sum(1 for g in gqs if g < thr_gq), len(gqs)),
            "mean_dp": mean(dps),
            "median_dp": percentile(dps, 0.5),
            "frac_dp_lt_thr": frac(sum(1 for d in dps if d < thr_dp), len(dps)),
            "frac_rnc_has_I": frac(sum(1 for e in subset if "I" in (e.rnc or "")), len(subset)),
            "label": label,
        }

    summary = {
        **meta,
        "switch_site_calls": {
            "n": qc.switch_n,
            "frac_low_gq": frac(qc.switch_low_gq, qc.switch_n),
            "frac_low_dp": frac(qc.switch_low_dp, qc.switch_n),
            "frac_rnc_I": frac(qc.switch_rnc_i, qc.switch_n),
            "mean_gq": mean(qc.switch_gq),
            "median_gq": percentile(qc.switch_gq, 0.5),
            "mean_dp": mean(qc.switch_dp),
            "median_dp": percentile(qc.switch_dp, 0.5),
        },
        "background_site_calls": {
            "n": qc.bg_n,
            "frac_low_gq": frac(qc.bg_low_gq, qc.bg_n),
            "frac_low_dp": frac(qc.bg_low_dp, qc.bg_n),
            "frac_rnc_I": frac(qc.bg_rnc_i, qc.bg_n),
            "mean_gq": mean(qc.bg_gq),
            "median_gq": percentile(qc.bg_gq, 0.5),
            "mean_dp": mean(qc.bg_dp),
            "median_dp": percentile(qc.bg_dp, 0.5),
        },
        "switches_all": qual_block("all", events),
        "switches_flicker": qual_block("flicker", flicker),
        "switches_sustained": qual_block("sustained", sustained),
        "by_gt": {},
    }

    def bin_block(bin_: GtBin) -> dict:
        return {
            "n": bin_.n,
            "frac_of_parent": None,
            "frac_low_gq": frac(bin_.low_gq, bin_.n),
            "frac_low_dp": frac(bin_.low_dp, bin_.n),
            "frac_sentinel": frac(bin_.sentinel, bin_.n),
            "frac_rnc_I": frac(bin_.rnc_i, bin_.n),
            "mean_gq": mean(bin_.gq),
            "median_gq": percentile(bin_.gq, 0.5),
            "mean_dp": mean(bin_.dp),
            "median_dp": percentile(bin_.dp, 0.5),
        }

    by_gt = {}
    for klass in GT_CLASSES:
        sw_bin = qc.switch_by_gt.get(klass, GtBin())
        bg_bin = qc.bg_by_gt.get(klass, GtBin())
        sw_d = bin_block(sw_bin)
        bg_d = bin_block(bg_bin)
        sw_d["frac_of_parent"] = frac(sw_bin.n, qc.switch_n)
        bg_d["frac_of_parent"] = frac(bg_bin.n, qc.bg_n)
        row = {
            "gt_class": klass,
            "switch": sw_d,
            "background": bg_d,
            "enrichment_frac_low_gq": None,
            "enrichment_frac_low_dp": None,
            "enrichment_frac_rnc_I": None,
        }
        for key in ("frac_low_gq", "frac_low_dp", "frac_rnc_I"):
            s = sw_d.get(key)
            b = bg_d.get(key)
            if s is not None and b is not None and b > 0:
                row[f"enrichment_{key}"] = s / b
        by_gt[klass] = row
    summary["by_gt"] = by_gt
    summary["switch_site_calls"]["frac_sentinel"] = frac(qc.switch_sentinel, qc.switch_n)
    summary["background_site_calls"]["frac_sentinel"] = frac(qc.bg_sentinel, qc.bg_n)
    s_sent = summary["switch_site_calls"]["frac_sentinel"]
    b_sent = summary["background_site_calls"]["frac_sentinel"]
    if s_sent is not None and b_sent is not None and b_sent > 0:
        summary["enrichment_frac_sentinel"] = s_sent / b_sent
    else:
        summary["enrichment_frac_sentinel"] = None

    # Enrichment: P(low GQ | switch) / P(low GQ | bg)
    sw = summary["switch_site_calls"]
    bg = summary["background_site_calls"]
    for key in ("frac_low_gq", "frac_low_dp", "frac_rnc_I"):
        s = sw.get(key)
        b = bg.get(key)
        if s is not None and b is not None and b > 0:
            summary[f"enrichment_{key}"] = s / b
        else:
            summary[f"enrichment_{key}"] = None

    def flanking_het_block(label: str, subset: list[SwitchEvent], *, gq_attr: str, dp_attr: str, gap_attr: str) -> dict:
        gqs = [getattr(e, gq_attr) for e in subset if getattr(e, gq_attr) is not None]
        dps = [getattr(e, dp_attr) for e in subset if getattr(e, dp_attr) is not None]
        gaps = [getattr(e, gap_attr) for e in subset if getattr(e, gap_attr) is not None]
        thr_gq = meta["gq_threshold"]
        thr_dp = meta["dp_threshold"]
        n_with = sum(1 for e in subset if getattr(e, gq_attr) is not None or getattr(e, dp_attr) is not None)
        return {
            "n": len(subset),
            "n_with_het": n_with,
            "frac_with_het": frac(n_with, len(subset)),
            "mean_gq": mean(gqs),
            "median_gq": percentile(gqs, 0.5),
            "frac_gq_lt_thr": frac(sum(1 for g in gqs if g < thr_gq), len(gqs)),
            "mean_dp": mean(dps),
            "median_dp": percentile(dps, 0.5),
            "frac_dp_lt_thr": frac(sum(1 for d in dps if d < thr_dp), len(dps)),
            "median_gap_bp": percentile(gaps, 0.5) if gaps else None,
            "label": label,
        }

    hom_ref_sw = [e for e in events if e.gt_class == "hom_ref"]
    het_sw = [e for e in events if e.gt_class == "het"]
    bg_het = qc.bg_by_gt.get("het", GtBin())
    nearest = {
        "all": flanking_het_block("nearest_het_all", events, gq_attr="nearest_het_gq", dp_attr="nearest_het_dp", gap_attr="nearest_het_gap"),
        "flicker": flanking_het_block("nearest_het_flicker", flicker, gq_attr="nearest_het_gq", dp_attr="nearest_het_dp", gap_attr="nearest_het_gap"),
        "sustained": flanking_het_block("nearest_het_sustained", sustained, gq_attr="nearest_het_gq", dp_attr="nearest_het_dp", gap_attr="nearest_het_gap"),
        "prev_het_all": flanking_het_block("prev_het_all", events, gq_attr="prev_het_gq", dp_attr="prev_het_dp", gap_attr="prev_het_gap"),
        "switch_at_hom_ref": flanking_het_block(
            "nearest_het_switch_at_hom_ref", hom_ref_sw,
            gq_attr="nearest_het_gq", dp_attr="nearest_het_dp", gap_attr="nearest_het_gap",
        ),
        "switch_at_het": flanking_het_block(
            "nearest_het_switch_at_het", het_sw,
            gq_attr="nearest_het_gq", dp_attr="nearest_het_dp", gap_attr="nearest_het_gap",
        ),
        "background_hets": {
            "n": bg_het.n,
            "median_gq": percentile(bg_het.gq, 0.5),
            "median_dp": percentile(bg_het.dp, 0.5),
            "frac_low_gq": frac(bg_het.low_gq, bg_het.n),
            "frac_low_dp": frac(bg_het.low_dp, bg_het.n),
        },
        "side_counts": {},
        "enrichment_frac_low_gq": None,
        "enrichment_frac_low_dp": None,
        "hom_ref_switch_enrichment_frac_low_gq": None,
        "hom_ref_switch_enrichment_frac_low_dp": None,
    }
    for side in ("current", "prev", "next", "none"):
        nearest["side_counts"][side] = sum(1 for e in events if e.nearest_het_side == side)
    bg_low_gq = nearest["background_hets"]["frac_low_gq"]
    bg_low_dp = nearest["background_hets"]["frac_low_dp"]
    for src_key, enr_gq, enr_dp in (
        ("all", "enrichment_frac_low_gq", "enrichment_frac_low_dp"),
        ("switch_at_hom_ref", "hom_ref_switch_enrichment_frac_low_gq", "hom_ref_switch_enrichment_frac_low_dp"),
    ):
        s_gq = nearest[src_key]["frac_gq_lt_thr"]
        s_dp = nearest[src_key]["frac_dp_lt_thr"]
        if s_gq is not None and bg_low_gq is not None and bg_low_gq > 0:
            nearest[enr_gq] = s_gq / bg_low_gq
        if s_dp is not None and bg_low_dp is not None and bg_low_dp > 0:
            nearest[enr_dp] = s_dp / bg_low_dp
    summary["nearest_het"] = nearest

    def rnc_block(xs: list[str]) -> dict:
        n = len(xs)
        codes = Counter(xs)
        return {
            "n": n,
            "frac_I": frac(sum(1 for x in xs if rnc_has_I(x)), n),
            "frac_clean": frac(sum(1 for x in xs if x in (".", "..")), n),
            "top_codes": dict(codes.most_common(12)),
        }

    hom_ref_rnc = [e.rnc for e in hom_ref_sw]
    summary["rnc"] = {
        "switch": rnc_block(qc.switch_rnc),
        "background": rnc_block(qc.bg_rnc),
        "switch_at_hom_ref": rnc_block(hom_ref_rnc),
        "switch_at_het": rnc_block([e.rnc for e in het_sw]),
        "nearest_het_all": rnc_block([e.nearest_het_rnc for e in events if e.nearest_het_side != "none"]),
        "nearest_het_switch_at_hom_ref": rnc_block(
            [e.nearest_het_rnc for e in hom_ref_sw if e.nearest_het_side != "none"]
        ),
    }
    sw_i = summary["rnc"]["switch"]["frac_I"]
    bg_i = summary["rnc"]["background"]["frac_I"]
    if sw_i is not None and bg_i is not None and bg_i > 0:
        summary["rnc"]["enrichment_frac_I"] = sw_i / bg_i
    else:
        summary["rnc"]["enrichment_frac_I"] = None
    hr_i = summary["rnc"]["switch_at_hom_ref"]["frac_I"]
    if hr_i is not None and bg_i is not None and bg_i > 0:
        summary["rnc"]["hom_ref_switch_enrichment_frac_I"] = hr_i / bg_i
    else:
        summary["rnc"]["hom_ref_switch_enrichment_frac_I"] = None
    if tracts is not None:
        summary["tracts"] = tract_stats(tracts, meta, events)
    else:
        summary["tracts"] = tract_stats([], meta, events)
    # Part 7.1: re-score switch rates after dropping poor GQ/DP/RNC calls.
    if meta.get("gq_dp_available", True) and events:
        summary["call_qc_filtered"] = call_qc_filtered_tract_metrics(events, meta)
    else:
        summary["call_qc_filtered"] = {
            "available": False,
            "reason": "GQ/DP unavailable or no switches",
        }
    return summary


def write_arrays_for_plots(
    out_dir: Path,
    qc: QcAccum,
    events: Optional[list[SwitchEvent]] = None,
    tracts: Optional[list[AncTract]] = None,
) -> None:
    """Small sidecar arrays the notebook can plot without re-scanning the VCF."""
    import csv

    def dump(name: str, xs: list[int]) -> None:
        p = out_dir / name
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["value"])
            # Cap for notebook friendliness
            for v in xs[:500_000]:
                w.writerow([v])

    dump("switch_gq.csv", qc.switch_gq)
    dump("switch_dp.csv", qc.switch_dp)
    dump("background_gq.csv", qc.bg_gq)
    dump("background_dp.csv", qc.bg_dp)

    def dump_str(name: str, xs: list[str]) -> None:
        p = out_dir / name
        with open(p, "w", newline="") as fh:
            w = csv.writer(fh)
            w.writerow(["value"])
            for v in xs[:500_000]:
                w.writerow([v])

    dump_str("switch_rnc.csv", qc.switch_rnc)
    dump_str("background_rnc.csv", qc.bg_rnc)

    by_gt = out_dir / "gq_dp_by_gt.csv"
    with open(by_gt, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["set", "gt_class", "metric", "value"])
        for set_name, bins in (("switch", qc.switch_by_gt), ("background", qc.bg_by_gt)):
            for klass, bin_ in bins.items():
                for metric, xs in (("gq", bin_.gq), ("dp", bin_.dp)):
                    for v in xs[:200_000]:
                        w.writerow([set_name, klass, metric, v])

    if tracts:
        dump("tract_length_bp.csv", [t.length_bp for t in tracts if t.length_bp > 0])
        dump(
            "complete_tract_length_bp.csv",
            [
                t.length_bp
                for t in tracts
                if t.length_bp > 0 and not t.left_censored and not t.right_censored
            ],
        )

    if events is None:
        return
    nearest_path = out_dir / "nearest_het_gq_dp.csv"
    with open(nearest_path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(
            [
                "switch_class",
                "gt_class",
                "side",
                "gap_bp",
                "gq",
                "dp",
                "prev_het_gq",
                "prev_het_dp",
                "prev_het_gap",
                "rnc",
                "nearest_het_rnc",
                "prev_het_rnc",
            ]
        )
        for e in events[:500_000]:
            w.writerow(
                [
                    "flicker" if e.is_flicker else "sustained",
                    e.gt_class,
                    e.nearest_het_side,
                    "" if e.nearest_het_gap is None else e.nearest_het_gap,
                    "" if e.nearest_het_gq is None else e.nearest_het_gq,
                    "" if e.nearest_het_dp is None else e.nearest_het_dp,
                    "" if e.prev_het_gq is None else e.prev_het_gq,
                    "" if e.prev_het_dp is None else e.prev_het_dp,
                    "" if e.prev_het_gap is None else e.prev_het_gap,
                    e.rnc,
                    e.nearest_het_rnc,
                    e.prev_het_rnc,
                ]
            )


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--vcf",
        type=Path,
        required=True,
        help="FLARE VCF/BCF with AN1/AN2 (GQ/DP/RNC optional)",
    )
    p.add_argument("--samples", type=Path, default=None, help="Optional sample ID list (one per line)")
    p.add_argument("--out-dir", type=Path, required=True)
    p.add_argument("--region", default=None, help="Optional chr:start-end passed to bcftools query -r")
    p.add_argument("--gq-threshold", type=int, default=20)
    p.add_argument("--dp-threshold", type=int, default=5)
    p.add_argument(
        "--flicker-max-bp",
        type=int,
        default=50_000,
        help="Max span for A→B→A flicker classification (default 50kb)",
    )
    p.add_argument(
        "--bg-keep-every",
        type=int,
        default=50,
        help="Keep every Nth non-switch sample-site for background QC (default 50)",
    )
    p.add_argument("--max-sites", type=int, default=None, help="Optional cap for smoke tests")
    p.add_argument(
        "--jobs",
        type=int,
        default=1,
        help="Parallel sample shards (each runs its own bcftools query). Default 1.",
    )
    p.add_argument(
        "--an-only",
        action="store_true",
        help="Query only AN1/AN2; skip GQ/DP/GT/RNC and het bookkeeping (fast rates/flicker).",
    )
    p.add_argument(
        "--summary-only",
        action="store_true",
        help="Skip full switches/tracts TSVs and plot arrays; write summary.json + switches_slim.tsv.gz.",
    )
    args = p.parse_args(argv)

    samples = load_sample_ids(args.samples, args.vcf)
    args.out_dir.mkdir(parents=True, exist_ok=True)
    jobs = max(1, int(args.jobs))

    print(
        f"Scanning {args.vcf} for {len(samples)} samples "
        f"(GQ<{args.gq_threshold}, DP<{args.dp_threshold}, flicker≤{args.flicker_max_bp} bp"
        + (f", region={args.region}" if args.region else "")
        + (f", jobs={jobs}" if jobs > 1 else "")
        + (", an-only" if args.an_only else "")
        + (", summary-only" if args.summary_only else "")
        + "; AN-only FLARE VCFs OK)",
        flush=True,
    )
    events, qc, meta, tracts = scan_vcf_parallel(
        args.vcf,
        samples,
        jobs=jobs,
        gq_threshold=args.gq_threshold,
        dp_threshold=args.dp_threshold,
        flicker_max_bp=args.flicker_max_bp,
        bg_keep_every=args.bg_keep_every,
        max_sites=args.max_sites,
        region=args.region,
        an_only=args.an_only,
        rates_only=args.summary_only or args.an_only,
    )
    events.sort(key=lambda e: (e.chrom, e.pos, e.sample, e.hap))

    slim_path = args.out_dir / "switches_slim.tsv.gz"
    write_switches_slim_tsv(slim_path, events)
    print(f"Wrote {slim_path} ({len(events)} switches)", flush=True)

    if not args.summary_only:
        switches_path = args.out_dir / "switches.tsv.gz"
        write_switches_tsv(switches_path, events)
        tracts_path = args.out_dir / "tracts.tsv.gz"
        write_tracts_tsv(tracts_path, tracts)
        print(f"Wrote {switches_path} ({len(events)} switches)", flush=True)
        print(f"Wrote {tracts_path} ({len(tracts)} tracts)", flush=True)
        write_arrays_for_plots(args.out_dir, qc, events, tracts)

    summary = summarize(qc, meta, events, tracts)
    summary_path = args.out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    print(f"Wrote {summary_path}", flush=True)
    tr = summary.get("tracts") or {}
    print(
        f"Tracts: span_mb={tr.get('span_mb')} "
        f"switches/hap/Mb={tr.get('switches_per_hap_per_mb')} "
        f"flicker_frac={tr.get('flicker_frac')} "
        f"flicker/hap/Mb={tr.get('flicker_per_hap_per_mb')} "
        f"implied_T={tr.get('implied_t_gen')} "
        f"median_len={tr.get('median_length_bp')}",
        flush=True,
    )
    cqf = summary.get("call_qc_filtered") or {}
    if cqf.get("n_switches_all") is not None:
        print(
            f"Call-QC filtered: fail_frac={cqf.get('frac_switches_fail_call_qc')} "
            f"switches/hap/Mb={cqf.get('switches_per_hap_per_mb')} "
            f"implied_T={cqf.get('implied_t_gen')} "
            f"(GQ>={cqf.get('gq_threshold')}, DP>={cqf.get('dp_threshold')}, drop_RNC_I)",
            flush=True,
        )
    sw = summary["switch_site_calls"]
    bg = summary["background_site_calls"]
    if summary.get("gq_dp_available"):
        print(
            f"Switch sites: n={sw['n']} median_GQ={sw['median_gq']} "
            f"frac_GQ<{args.gq_threshold}={sw['frac_low_gq']}",
            flush=True,
        )
        print(
            f"Background:   n={bg['n']} median_GQ={bg['median_gq']} "
            f"frac_GQ<{args.gq_threshold}={bg['frac_low_gq']}",
            flush=True,
        )
        print(
            f"Enrichment low-GQ={summary.get('enrichment_frac_low_gq')} "
            f"low-DP={summary.get('enrichment_frac_low_dp')} "
            f"RNC-I={summary.get('enrichment_frac_rnc_I')} "
            f"GLnexus-homref-sentinel={summary.get('enrichment_frac_sentinel')}",
            flush=True,
        )
        het = summary.get("by_gt", {}).get("het", {})
        het_sw = het.get("switch", {})
        print(
            f"Het only: switch n={het_sw.get('n')} median_GQ={het_sw.get('median_gq')} "
            f"median_DP={het_sw.get('median_dp')} "
            f"enrichment_low_GQ={het.get('enrichment_frac_low_gq')}",
            flush=True,
        )
        nh = summary.get("nearest_het", {})
        nh_all = nh.get("all", {})
        nh_hr = nh.get("switch_at_hom_ref", {})
        print(
            f"Nearest het: n_with={nh_all.get('n_with_het')} median_GQ={nh_all.get('median_gq')} "
            f"median_gap={nh_all.get('median_gap_bp')} "
            f"enrichment_low_GQ={nh.get('enrichment_frac_low_gq')}",
            flush=True,
        )
        print(
            f"Nearest het (switch at hom_ref): n={nh_hr.get('n')} "
            f"median_GQ={nh_hr.get('median_gq')} median_DP={nh_hr.get('median_dp')} "
            f"median_gap={nh_hr.get('median_gap_bp')} "
            f"enrichment_low_GQ={nh.get('hom_ref_switch_enrichment_frac_low_gq')}",
            flush=True,
        )
        rnc = summary.get("rnc", {})
        print(
            f"RNC-I: switch={rnc.get('switch', {}).get('frac_I')} "
            f"background={rnc.get('background', {}).get('frac_I')} "
            f"enrichment={rnc.get('enrichment_frac_I')} "
            f"hom_ref_switch={rnc.get('switch_at_hom_ref', {}).get('frac_I')}",
            flush=True,
        )
    else:
        print(
            "GQ/DP unavailable (gq_dp_available=false); skipping quality summary lines.",
            flush=True,
        )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
