#!/usr/bin/env python3
"""Put TRGT alleles on the phased-callset haplotypes, by a read-level join.

TRGT genotypes are unphased: ``GT`` is ``/``-separated and the order of the
per-haplotype FORMAT fields (``AL``, ``ALLR``, ``SD``, ``MC``, ``MS``, ``AP``,
``AM``) is arbitrary. A per-ancestry repeat-dosage test needs each allele on
the *same* haplotype index that carries the FLARE ancestry, which is the
haplotype of the joint phased callset (HiPhase -> SHAPEIT4), not TRGT's.

This script does not infer orientation from any length estimate. It joins two
read-level tables that both key on read name::

    read -> haplotype   whatshap haplotag --output-haplotag-list, run against
                        the phased callset, so haplotype 1/2 are *our* indices
    read -> slot        TRGT's spanning-read BAM, which has carried
                        read-to-allele assignments since TRGT 0.5.0. The tag
                        holds the GENOTYPE SLOT (0 or 1), not the ALT index
                        and not an allele length: at a 1/2 locus the slots
                        carry alleles 1 and 2 while reads still say 0 and 1.

A single HiFi molecule physically carries both the repeat allele and the
flanking heterozygous SNVs that decide its haplotype, so the join is a
measurement rather than a model. The only thing it assumes beyond the rest of
the analysis is that those two tools agree about read names.

Orientation can only matter where the two genotype slots hold different allele
indices. ``0/0``, ``1/1`` and any call whose slots are equal are emitted
unchanged with ``phase_source = homozygous``: swapping identical alleles is a
no-op, not a resolved phase.

The allele-assignment tag in the spanning BAM is **not** hardcoded, because it
is version-dependent. Run with ``--census-only`` first: that writes the tag
census and exits, so you can see which tags the file actually carries before
committing to one. ``--allele-tag`` is required otherwise, and a tag that is
absent is a hard error with the census attached, never a silent empty join.

Outputs
-------
``--out-alleles``  one row per locus per haplotype: locus_id, chrom, start,
                   end, hap (1 or 2, our frame), al_bp, change_bp,
                   change_units, n_reads, purity, phase_source
``--out-vcf``      the input records with the genotype slots reordered so slot
                   0 is haplotype 1, ``GT`` pipe-separated, and INFO
                   ``TRGT_PHASE`` recording how each record was resolved
``--out-summary``  counts, the tag census, and the per-locus resolution tally
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import json
import re
import subprocess
import sys
from collections import Counter, defaultdict
from pathlib import Path
from typing import Iterator, Optional, TextIO

# FORMAT keys whose value is one comma-separated entry per genotype slot. These
# must be reordered with GT or the record silently misreports which haplotype
# carries which length. MC and MS use '_' within a slot and ',' between slots.
PER_SLOT_KEYS = ("AL", "ALLR", "SD", "MC", "MS", "AP", "AM")
GT_SPLIT = re.compile(r"[/|]")
TAG_RE = re.compile(r"^([A-Za-z][A-Za-z0-9]):([AifZHB]):(.*)$")
HAP_RE = re.compile(r"^(?:H|HP|hap)?([12])$", re.IGNORECASE)


def open_text(path: str, mode: str = "rt") -> TextIO:
    if str(path).endswith(".gz"):
        return gzip.open(path, mode)
    return open(path, mode)


def read_haplotag_list(path: str) -> tuple[dict[str, int], Counter]:
    """read name -> haplotype index (1 or 2), from whatshap's haplotag list.

    whatshap writes ``#readname haplotype phaseset chromosome`` with haplotype
    ``H1`` / ``H2`` / ``none``. HiPhase's own ``--haplotag-file`` uses a
    ``haplotag`` column holding ``1`` / ``2``; both are accepted, by name where
    a header exists and by position otherwise. Untagged reads are counted and
    dropped: a read with no haplotype carries no orientation information.
    """
    tally: Counter = Counter()
    out: dict[str, int] = {}
    with open_text(path) as fh:
        header: Optional[list[str]] = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            fields = line.split("\t")
            if header is None and (line.startswith("#") or _looks_like_header(fields)):
                header = [f.lstrip("#").strip() for f in fields]
                continue
            if header:
                row = dict(zip(header, fields))
                name = row.get("readname") or row.get("read_name") or fields[0]
                hap_raw = row.get("haplotype") or row.get("haplotag") or ""
            else:
                name, hap_raw = fields[0], (fields[1] if len(fields) > 1 else "")
            m = HAP_RE.match(hap_raw.strip())
            if not m:
                tally["read_untagged"] += 1
                continue
            out[name] = int(m.group(1))
            tally[f"read_hap{m.group(1)}"] += 1
    return out, tally


def _looks_like_header(fields: list[str]) -> bool:
    lowered = {f.lstrip("#").strip().lower() for f in fields}
    return bool(lowered & {"readname", "read_name", "haplotype", "haplotag"})


def vcf_lines(path: str) -> Iterator[str]:
    """Stream a VCF. bcftools handles BCF and gs:// URIs; plain text is read directly."""
    if path.endswith((".bcf", ".bcf.gz")) or path.startswith("gs://"):
        proc = subprocess.Popen(["bcftools", "view", path], stdout=subprocess.PIPE, text=True,
                                bufsize=1 << 20)
        assert proc.stdout is not None
        try:
            yield from proc.stdout
        finally:
            proc.stdout.close()
            if proc.wait() != 0:
                raise SystemExit(f"bcftools view failed on {path}")
    else:
        with open_text(path) as fh:
            yield from fh


def info_get(info: str, key: str) -> Optional[str]:
    for field in info.split(";"):
        if field == key:
            return ""
        if field.startswith(key + "="):
            return field[len(key) + 1:]
    return None


class LocusIndex:
    """Catalog intervals from the TRGT VCF, for assigning a read to a locus.

    The spanning BAM holds read portions placed at the repeat, so overlap with
    the record's [POS, END] identifies the locus. Where intervals overlap, the
    locus with the largest overlap wins and the read is counted as ambiguous.
    """

    def __init__(self) -> None:
        self.by_chrom: dict[str, list[tuple[int, int, str]]] = defaultdict(list)
        self._starts: dict[str, list[int]] = {}

    def add(self, chrom: str, start: int, end: int, trid: str) -> None:
        self.by_chrom[chrom].append((start, end, trid))

    def finalize(self) -> None:
        for chrom, items in self.by_chrom.items():
            items.sort()
            self._starts[chrom] = [s for s, _, _ in items]

    def lookup(self, chrom: str, start: int, end: int) -> tuple[Optional[str], bool]:
        items = self.by_chrom.get(chrom)
        if not items:
            return None, False
        starts = self._starts[chrom]
        i = bisect.bisect_right(starts, end)
        best, best_ovl, hits = None, 0, 0
        for s, e, trid in reversed(items[:i]):
            if e < start:
                # Intervals are start-sorted; a long locus further back can still
                # reach, so scan a bounded window rather than breaking at once.
                if s < start - 100_000:
                    break
                continue
            ovl = min(end, e) - max(start, s)
            if ovl <= 0:
                continue
            hits += 1
            if ovl > best_ovl:
                best, best_ovl = trid, ovl
        return best, hits > 1


def tag_census(bam: str, limit: int, region: str = "") -> tuple[Counter, dict[str, int]]:
    """Count how often each SAM tag appears, over the first ``limit`` records."""
    tags: Counter = Counter()
    n = 0
    for name, chrom, pos, ln, rec_tags in iter_bam(bam, region):
        n += 1
        for tag in rec_tags:
            tags[tag] += 1
        if limit and n >= limit:
            break
    return tags, {"records_scanned": n}


def iter_bam(bam: str, region: str = "") -> Iterator[tuple[str, str, int, int, dict[str, str]]]:
    """(read name, chrom, start, end, tags) per alignment, via samtools view.

    stderr is captured and surfaced: a missing index, a permissions failure and
    an empty region otherwise all look like zero reads.

    A caller that stops early (the tag census takes only the first N records)
    closes this generator, samtools then takes SIGPIPE, and its non-zero exit
    says nothing about whether the records read were sound. So the exit status
    is only checked when the stream was consumed to the end.
    """
    cmd = ["samtools", "view", bam] + ([region] if region else [])
    proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            bufsize=1 << 20)
    assert proc.stdout is not None and proc.stderr is not None
    completed = False
    try:
        for line in proc.stdout:
            f = line.rstrip("\n").split("\t")
            if len(f) < 11:
                continue
            tags: dict[str, str] = {}
            for field in f[11:]:
                m = TAG_RE.match(field)
                if m:
                    tags[m.group(1)] = m.group(3)
            start = int(f[3]) - 1
            yield f[0], f[2], start, start + len(f[9]), tags
        completed = True
    finally:
        if not completed:
            proc.kill()
        proc.stdout.close()
        err = proc.stderr.read()
        proc.stderr.close()
        rc = proc.wait()
        if completed and rc != 0:
            raise SystemExit(f"samtools view failed on {bam}:\n{err.strip()}")


def collect_votes(bam: str, loci: LocusIndex, haps: dict[str, int], allele_tag: str,
                  region: str, trid_tag: str = "", known: frozenset = frozenset()
                  ) -> tuple[dict[str, Counter], Counter]:
    """(locus, allele index, haplotype) vote counts from the spanning reads."""
    votes: dict[str, Counter] = defaultdict(Counter)
    tally: Counter = Counter()
    seen_tag = False
    for name, chrom, start, end, tags in iter_bam(bam, region):
        tally["spanning_reads"] += 1
        raw = tags.get(allele_tag)
        if raw is None:
            tally["read_no_allele_tag"] += 1
            continue
        seen_tag = True
        try:
            slot = int(raw)
        except ValueError:
            tally["read_allele_tag_unparsed"] += 1
            continue
        if slot not in (0, 1):
            # TRGT writes AL:i:0 / AL:i:1 -- the genotype SLOT the read was
            # assigned to, not the ALT index and not a length. Anything else
            # means the wrong tag was named.
            tally["read_slot_out_of_range"] += 1
            continue
        hap = haps.get(name)
        if hap is None:
            tally["read_not_haplotagged"] += 1
            continue
        # TRGT stamps the locus id on every spanning record, so prefer it. The
        # catalog has overlapping definitions (chr20:10,002,372-10,002,525 and
        # 10,002,512-10,002,565 among others), which makes assignment by
        # interval overlap genuinely ambiguous; the tag is exact.
        trid = tags.get(trid_tag) if trid_tag else None
        if trid is not None:
            trid = trid.split(",")[0]
            if trid not in known:
                tally["read_locus_not_in_vcf"] += 1
                continue
            tally["read_locus_from_tag"] += 1
        else:
            if trid_tag:
                tally["read_no_locus_tag"] += 1
            trid, ambiguous = loci.lookup(chrom, start, end)
            if trid is None:
                tally["read_no_locus"] += 1
                continue
            if ambiguous:
                tally["read_ambiguous_locus"] += 1
        votes[trid][(slot, hap)] += 1
        tally["read_joined"] += 1
    if not seen_tag and tally["spanning_reads"]:
        raise SystemExit(
            f"error: no alignment in {bam} carries tag {allele_tag!r}. "
            "Re-run with --census-only to list the tags this file actually has."
        )
    return votes, tally


def assign(slots: list[str], votes: Counter, min_reads: int, min_purity: float
           ) -> tuple[Optional[list[int]], str, dict[int, tuple[int, float]]]:
    """Order the genotype slots so slot 0 is haplotype 1.

    ``slots`` holds the allele index per genotype slot, as written by TRGT;
    the votes are keyed by slot POSITION, which is what the read tag carries.
    Returns the permutation to apply, why, and per-slot (read count, purity).
    Homozygous and single-slot calls need no evidence: the permutation is the
    identity and nothing is claimed about phase.
    """
    support: dict[int, tuple[int, float]] = {}
    if len(slots) < 2 or slots[0] == slots[1]:
        return list(range(len(slots))), "homozygous", support
    hap_of: dict[int, int] = {}
    for pos in range(len(slots)):
        counts = Counter({h: votes.get((pos, h), 0) for h in (1, 2)})
        total = sum(counts.values())
        hap, top = counts.most_common(1)[0]
        purity = (top / total) if total else 0.0
        support[pos] = (total, purity)
        if total < min_reads or purity < min_purity or top == 0:
            return None, "insufficient_reads", support
        hap_of[pos] = hap
    if hap_of[0] == hap_of[1]:
        return None, "conflict", support
    perm = sorted(range(len(slots)), key=lambda p: hap_of[p])
    return perm, "reads", support


def reorder_sample(fmt: str, sample: str, perm: list[int], phased: bool) -> str:
    """Apply ``perm`` to GT and every per-haplotype FORMAT field."""
    keys = fmt.split(":")
    values = sample.split(":")
    values += ["."] * (len(keys) - len(values))
    out: list[str] = []
    for key, value in zip(keys, values):
        if key == "GT":
            alleles = GT_SPLIT.split(value)
            if len(alleles) == len(perm):
                sep = "|" if phased else "/"
                value = sep.join(alleles[p] for p in perm)
        elif key in PER_SLOT_KEYS and value not in (".", ""):
            parts = value.split(",")
            if len(parts) == len(perm):
                value = ",".join(parts[p] for p in perm)
        out.append(value)
    return ":".join(out)


def units(change_bp: int, period: int) -> float:
    return change_bp / period if period else float(change_bp)


def run(args: argparse.Namespace) -> dict:
    census, scan = tag_census(args.spanning_bam, args.census_records, args.region)
    summary: dict = {
        "trgt_vcf": args.trgt_vcf,
        "spanning_bam": args.spanning_bam,
        "haplotag_list": args.haplotag_list,
        "allele_tag": args.allele_tag,
        "region": args.region or None,
        "min_reads": args.min_reads,
        "min_purity": args.min_purity,
        "tag_census": dict(sorted(census.items())),
        "tag_census_records": scan["records_scanned"],
    }
    if args.census_only:
        summary["status"] = "census_only"
        return summary

    haps, hap_tally = read_haplotag_list(args.haplotag_list)
    summary["haplotag"] = dict(hap_tally)
    summary["haplotag_reads"] = len(haps)

    records: list[tuple[list[str], str, int, int, str, int]] = []
    header: list[str] = []
    loci = LocusIndex()
    for line in vcf_lines(args.trgt_vcf):
        if line.startswith("#"):
            header.append(line.rstrip("\n"))
            continue
        f = line.rstrip("\n").split("\t")
        if len(f) < 10:
            continue
        chrom, pos, info = f[0], int(f[1]), f[7]
        trid = info_get(info, "TRID") or f"{chrom}:{pos}"
        trid = trid.split(",")[0]
        end_raw = info_get(info, "END")
        end = int(end_raw) if end_raw else pos
        motifs = (info_get(info, "MOTIFS") or "").split(",")
        period = len(motifs[0]) if motifs and motifs[0] else 1
        loci.add(chrom, pos - 1, end, trid)
        records.append((f, trid, pos, end, motifs[0] if motifs else "", period))
    loci.finalize()
    summary["loci"] = len(records)

    votes, read_tally = collect_votes(args.spanning_bam, loci, haps, args.allele_tag,
                                      args.region, args.trid_tag,
                                      frozenset(r[1] for r in records))
    summary["reads"] = dict(read_tally)

    resolution: Counter = Counter()
    alleles_out = open_text(args.out_alleles, "wt") if args.out_alleles else None
    vcf_out = open_text(args.out_vcf, "wt") if args.out_vcf else None
    if alleles_out:
        alleles_out.write("\t".join((
            "locus_id", "chrom", "start", "end", "period", "motif", "hap",
            "al_bp", "change_bp", "change_units", "n_reads", "purity", "phase_source")) + "\n")
    if vcf_out:
        for line in header:
            if line.startswith("##fileformat"):
                vcf_out.write(line + "\n")
                vcf_out.write('##INFO=<ID=TRGT_PHASE,Number=1,Type=String,Description='
                              '"How the genotype slots were ordered onto the phased-callset '
                              'haplotypes: reads, homozygous, conflict, insufficient_reads">\n')
            else:
                vcf_out.write(line + "\n")

    for f, trid, pos, end, motif, period in records:
        fmt, sample = f[8], f[9]
        keys = fmt.split(":")
        values = sample.split(":")
        gt = values[keys.index("GT")] if "GT" in keys else "./."
        al = values[keys.index("AL")] if "AL" in keys else ""
        slots_raw = GT_SPLIT.split(gt)
        slots = [int(a) for a in slots_raw if a.isdigit()]
        if len(slots) != len(slots_raw):
            resolution["missing_gt"] += 1
            perm, why, support = None, "missing_gt", {}
        else:
            perm, why, support = assign(slots, votes.get(trid, Counter()),
                                        args.min_reads, args.min_purity)
            resolution[why] += 1
        phased = why == "reads"
        if vcf_out:
            use = perm if perm is not None else list(range(len(slots_raw)))
            f = list(f)
            f[7] = f[7] + f";TRGT_PHASE={why}"
            f[9] = reorder_sample(fmt, sample, use, phased)
            vcf_out.write("\t".join(f) + "\n")
        if alleles_out and perm is not None and why in ("reads", "homozygous"):
            ref_len = end - pos + 1
            lengths = al.split(",") if al else []
            for hap_idx, slot in enumerate(perm, start=1):
                if slot >= len(lengths) or not lengths[slot].strip().isdigit():
                    continue
                al_bp = int(lengths[slot])
                change = al_bp - ref_len
                n_reads, purity = support.get(slot, (0, 1.0 if why == "homozygous" else 0.0))
                alleles_out.write("\t".join((
                    trid, f[0], str(pos - 1), str(end), str(period), motif, str(hap_idx),
                    str(al_bp), str(change), f"{units(change, period):.6g}",
                    str(n_reads), f"{purity:.4g}", why)) + "\n")
    if alleles_out:
        alleles_out.close()
    if vcf_out:
        vcf_out.close()

    summary["resolution"] = dict(sorted(resolution.items()))
    n_orientable = resolution.get("reads", 0) + resolution.get("conflict", 0) \
        + resolution.get("insufficient_reads", 0)
    summary["n_orientation_relevant"] = n_orientable
    summary["frac_resolved_of_relevant"] = (
        resolution.get("reads", 0) / n_orientable if n_orientable else None
    )
    summary["status"] = "ok"
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--trgt-vcf", required=True, help="TRGT VCF for one sample")
    p.add_argument("--spanning-bam", required=True, help="TRGT spanning-read BAM for the same run")
    p.add_argument("--haplotag-list", help="whatshap haplotag --output-haplotag-list TSV")
    p.add_argument("--allele-tag", help="SAM tag holding the read's allele index "
                                        "(see --census-only; not defaulted, it is version-dependent)")
    p.add_argument("--trid-tag", default="TR",
                   help="SAM tag holding the locus id (TRGT writes TR on every spanning "
                        "record). Pass '' to assign reads to loci by interval overlap "
                        "instead, which is ambiguous where catalog definitions overlap.")
    p.add_argument("--region", default="", help="chrom:start-end; restricts the BAM read")
    p.add_argument("--min-reads", type=int, default=2,
                   help="Haplotagged spanning reads required per allele to claim phase")
    p.add_argument("--min-purity", type=float, default=0.8,
                   help="Fraction of an allele's reads that must agree on one haplotype")
    p.add_argument("--census-records", type=int, default=0,
                   help="Alignments to scan for the tag census (0 = all, the default: a "
                        "truncated census misreports tag coverage, since the tags are not "
                        "uniformly distributed along the contig)")
    p.add_argument("--census-only", action="store_true",
                   help="Write the tag census and stop, without joining")
    p.add_argument("--out-alleles", default=None)
    p.add_argument("--out-vcf", default=None)
    p.add_argument("--out-summary", required=True)
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    if not args.census_only:
        missing = [n for n, v in (("--haplotag-list", args.haplotag_list),
                                  ("--allele-tag", args.allele_tag)) if not v]
        if missing:
            raise SystemExit(f"error: {', '.join(missing)} required unless --census-only")
    summary = run(args)
    Path(args.out_summary).write_text(json.dumps(summary, indent=2, sort_keys=False) + "\n")
    json.dump(summary, sys.stdout, indent=2)
    sys.stdout.write("\n")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
