#!/usr/bin/env python3
"""Detect FLARE ancestry switches and compare GQ / DP at those sites.

Streams an annotated VCF/BCF from ``propagate-annotations`` (FORMAT includes
AN1, AN2, GQ, DP). Inputs must be sorted by chromosome position.

Example::

    python3 scripts/flare_switch_qc.py \\
      --vcf annotated.chr22.bcf \\
      --samples analysis_samples.txt \\
      --out-dir flare_switch_qc_chr22 \\
      --gq-threshold 20 \\
      --dp-threshold 5
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import os
import shutil
import subprocess
import sys
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable, Iterator, Optional, TextIO


ANCESTRY = {0: "eas", 1: "amr", 2: "eur", 3: "afr", 4: "sas"}


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


@dataclass
class SampleTracker:
    hap1: HapState = field(default_factory=HapState)
    hap2: HapState = field(default_factory=HapState)
    switches: list[SwitchEvent] = field(default_factory=list)


@dataclass
class QcAccum:
    """On-the-fly histograms for switch vs background quality."""

    switch_gq: list[int] = field(default_factory=list)
    switch_dp: list[int] = field(default_factory=list)
    bg_gq: list[int] = field(default_factory=list)
    bg_dp: list[int] = field(default_factory=list)
    switch_low_gq: int = 0
    switch_low_dp: int = 0
    switch_rnc_i: int = 0
    switch_n: int = 0
    bg_low_gq: int = 0
    bg_low_dp: int = 0
    bg_rnc_i: int = 0
    bg_n: int = 0
    sites_scanned: int = 0
    sample_sites: int = 0


def parse_int(tok: str) -> Optional[int]:
    if tok in ("", ".", "./."):
        return None
    try:
        return int(tok)
    except ValueError:
        return None


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


def iter_annotated_sites(
    vcf: Path,
    samples: list[str],
) -> Iterator[tuple[str, int, str, str, list[tuple[Optional[int], Optional[int], Optional[int], Optional[int], str]]]]:
    """Yield (chrom, pos, ref, alt, per_sample[(gq, dp, an1, an2, rnc), ...])."""
    fmt = r"%CHROM\t%POS\t%REF\t%ALT[\t%GQ\t%DP\t%AN1\t%AN2]\n"
    fields_per = 4
    cmd = [
        bcftools_bin(),
        "query",
        "-f",
        fmt,
        "-s",
        ",".join(samples),
        str(vcf),
    ]
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
            expect = 4 + fields_per * n_samp
            if len(parts) < expect:
                raise RuntimeError(
                    f"line {line_no}: expected {expect} fields, got {len(parts)}. "
                    "Is this a propagate-annotations VCF with GQ/DP/AN1/AN2?"
                )
            chrom, pos_s, ref, alt = parts[0], parts[1], parts[2], parts[3]
            pos = int(pos_s)
            per = []
            base = 4
            for _ in range(n_samp):
                gq = parse_int(parts[base])
                dp = parse_int(parts[base + 1])
                an1 = parse_int(parts[base + 2])
                an2 = parse_int(parts[base + 3])
                per.append((gq, dp, an1, an2, "."))
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
) -> Optional[SwitchEvent]:
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
            ref=ref,
            alt=alt,
            anc_from=state.anc,
            anc_to=anc,
            gq=gq,
            dp=dp,
            rnc=rnc,
            prev_pos=state.pos,
            bp_gap=pos - state.pos,
        )
        tracker.switches.append(ev)
    state.anc = anc
    state.pos = pos
    state.gq = gq
    state.dp = dp
    return ev


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
) -> None:
    def low_gq(v: Optional[int]) -> bool:
        return v is not None and v < gq_thr

    def low_dp(v: Optional[int]) -> bool:
        return v is not None and v < dp_thr

    rnc_bad = "I" in (rnc or "")
    if is_switch:
        qc.switch_n += 1
        if gq is not None:
            qc.switch_gq.append(gq)
        if dp is not None:
            qc.switch_dp.append(dp)
        if low_gq(gq):
            qc.switch_low_gq += 1
        if low_dp(dp):
            qc.switch_low_dp += 1
        if rnc_bad:
            qc.switch_rnc_i += 1
        return

    bg_counter[0] += 1
    if bg_counter[0] % bg_keep_every != 0:
        return
    qc.bg_n += 1
    if gq is not None:
        if len(qc.bg_gq) < 2_000_000:
            qc.bg_gq.append(gq)
    if dp is not None:
        if len(qc.bg_dp) < 2_000_000:
            qc.bg_dp.append(dp)
    if low_gq(gq):
        qc.bg_low_gq += 1
    if low_dp(dp):
        qc.bg_low_dp += 1
    if rnc_bad:
        qc.bg_rnc_i += 1


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


def scan_vcf(
    vcf: Path,
    samples: list[str],
    *,
    gq_threshold: int,
    dp_threshold: int,
    flicker_max_bp: int,
    bg_keep_every: int,
    max_sites: Optional[int] = None,
) -> tuple[list[SwitchEvent], QcAccum, dict]:
    trackers = {s: SampleTracker() for s in samples}
    qc = QcAccum()
    bg_counter = [0]
    chrom_sites: dict[str, int] = defaultdict(int)

    for chrom, pos, ref, alt, per in iter_annotated_sites(vcf, samples):
        qc.sites_scanned += 1
        chrom_sites[chrom] += 1
        for sample, (gq, dp, an1, an2, rnc) in zip(samples, per):
            tr = trackers[sample]
            qc.sample_sites += 1
            ev1 = maybe_record_switch(tr, sample, 1, chrom, pos, ref, alt, an1, gq, dp, rnc)
            ev2 = maybe_record_switch(tr, sample, 2, chrom, pos, ref, alt, an2, gq, dp, rnc)
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
            )
        if max_sites is not None and qc.sites_scanned >= max_sites:
            break

    all_switches: list[SwitchEvent] = []
    for tr in trackers.values():
        mark_flickers(tr.switches, flicker_max_bp)
        all_switches.extend(tr.switches)

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
    }
    return all_switches, qc, meta


def write_switches_tsv(path: Path, events: list[SwitchEvent]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with open_text(path, "wt") as fh:
        fh.write(
            "sample\thap\tchrom\tpos\tref\talt\tanc_from\tanc_to\t"
            "anc_from_label\tanc_to_label\tgq\tdp\trnc\tprev_pos\tbp_gap\t"
            "is_flicker\tflicker_span_bp\n"
        )
        for e in events:
            fh.write(
                f"{e.sample}\t{e.hap}\t{e.chrom}\t{e.pos}\t{e.ref}\t{e.alt}\t"
                f"{e.anc_from}\t{e.anc_to}\t"
                f"{ANCESTRY.get(e.anc_from, '?')}\t{ANCESTRY.get(e.anc_to, '?')}\t"
                f"{'' if e.gq is None else e.gq}\t"
                f"{'' if e.dp is None else e.dp}\t"
                f"{e.rnc}\t{e.prev_pos}\t{e.bp_gap}\t"
                f"{int(e.is_flicker)}\t"
                f"{'' if e.flicker_span_bp is None else e.flicker_span_bp}\n"
            )


def summarize(qc: QcAccum, meta: dict, events: list[SwitchEvent]) -> dict:
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
    }

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
    return summary


def write_arrays_for_plots(out_dir: Path, qc: QcAccum) -> None:
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


def main(argv: Optional[list[str]] = None) -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--vcf", type=Path, required=True, help="Annotated VCF/VCF.gz from propagate-annotations")
    p.add_argument("--samples", type=Path, default=None, help="Optional sample ID list (one per line)")
    p.add_argument("--out-dir", type=Path, required=True)
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
    args = p.parse_args(argv)

    samples = load_sample_ids(args.samples, args.vcf)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    print(
        f"Scanning {args.vcf} for {len(samples)} samples "
        f"(GQ<{args.gq_threshold}, DP<{args.dp_threshold}, flicker≤{args.flicker_max_bp} bp)",
        flush=True,
    )
    events, qc, meta = scan_vcf(
        args.vcf,
        samples,
        gq_threshold=args.gq_threshold,
        dp_threshold=args.dp_threshold,
        flicker_max_bp=args.flicker_max_bp,
        bg_keep_every=args.bg_keep_every,
        max_sites=args.max_sites,
    )
    events.sort(key=lambda e: (e.chrom, e.pos, e.sample, e.hap))

    switches_path = args.out_dir / "switches.tsv.gz"
    write_switches_tsv(switches_path, events)
    summary = summarize(qc, meta, events)
    summary_path = args.out_dir / "summary.json"
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")
    write_arrays_for_plots(args.out_dir, qc)

    print(f"Wrote {switches_path} ({len(events)} switches)", flush=True)
    print(f"Wrote {summary_path}", flush=True)
    sw = summary["switch_site_calls"]
    bg = summary["background_site_calls"]
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
        f"RNC-I={summary.get('enrichment_frac_rnc_I')}",
        flush=True,
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
