#!/usr/bin/env python3
"""
One repeat-unit dosage per haplotype per catalog locus, from an integrated callset.

At a tandem repeat, the integrated SNV/indel/SV callset splits one locus over
several records (GLnexus indels, merged SV records). Every length-changing
record that overlaps a catalog locus (``--catalog-bed``, e.g. the TRExplorer
intervals) is assigned to that locus, and each haplotype's dosage is the sum
of the signed length changes it carries::

    dosage_bp(h) = sum over records r in the locus of len(ALT_r[a_h]) - len(REF_r)
    dosage_units(h) = dosage_bp(h) / period          (period = locus motif length)

Signs come from REF/ALT lengths (symbolic ALTs: SVTYPE, DEL < 0), so
sequence-resolved deletions, mixed INS/DEL records and impure alleles are all
kept; RU_TEST is not required.

A haplotype is missing at a locus when any of its records is missing, and a
sample is missing when two or more of its records at the locus are unphased
heterozygous length changes (their haplotype assignment is ambiguous). One
unphased het plus any number of homozygous records is resolved: haplotype
order is arbitrary, so the het change goes on haplotype 1.

Small-variant and SV callsets both report some events. With ``--small-vcf``
and ``--sv-vcf``, an allele counts only from the callset that owns its size
(``|change| < --split-bp`` small, otherwise SV); the other copy is REF.

Outputs (``--out-alleles`` and ``--out-summary`` hold no sample ids):

  --out-alleles  locus_id chrom start end period motif n_motifs n_records
                 n_hap n_hap_missing dosage_bp dosage_units n_hap_allele
                 (one row per distinct dosage, REF-length allele included)
  --out-vcf      one record per polymorphic locus: ALT ``<TR:+Nbp>``, INFO
                 RU_TEST / RU_DOSAGE (signed units per ALT) / SVLEN / PERIOD /
                 MOTIF, FORMAT GT (and AN1:AN2 with --with-ancestry), for
                 ``write_admixed_dosage_vcf.py``
  --out-summary  record / locus / haplotype tallies, including haplotypes that
                 carry a <50 bp and a >=50 bp record of the same sign at one
                 locus (possible double counts between the indel and SV calls)

Example::

    python3 aggregate_repeat_loci.py \\
      --small-vcf glnexus.chr22.vcf.gz --sv-vcf v3_main.bcf --apply-filters PASS,. \\
      --catalog-bed trexplorer.bed.gz --chrom chr22 --samples participants.txt \\
      --out-alleles chr22.alleles.tsv.gz --out-summary chr22.summary.json
"""

from __future__ import annotations

import argparse
import bisect
import gzip
import heapq
import json
import os
import subprocess
import sys
import tempfile
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterator, Optional, Sequence, TextIO

import numpy as np

SV_MIN_BP = 50
GAIN_TYPES = ("INS", "DUP", "CNV")


@dataclass
class Locus:
    chrom: str
    start: int  # 0-based
    end: int  # half-open
    lid: str
    motif: str
    n_motifs: int
    period_bp: int = 0

    @property
    def period(self) -> int:
        return self.period_bp or max(len(self.motif), 1)


@dataclass
class Record:
    pos: int
    start0: int
    end0: int
    rid: str
    alt_bp: np.ndarray  # signed bp change per allele index (0 = REF)
    a1: np.ndarray  # int16 allele index per sample, -1 missing
    a2: np.ndarray
    unphased_het: np.ndarray  # bool
    an1: Optional[np.ndarray] = None
    an2: Optional[np.ndarray] = None
    uncatalogued_period: Optional[int] = None
    uncatalogued_motif: Optional[str] = None


@dataclass
class LocusState:
    locus: Locus
    n: int
    with_ancestry: bool
    n_records: int = 0
    sum1: np.ndarray = field(init=False)
    sum2: np.ndarray = field(init=False)
    miss1: np.ndarray = field(init=False)
    miss2: np.ndarray = field(init=False)
    nonref_records: np.ndarray = field(init=False)
    unphased_hets: np.ndarray = field(init=False)
    carried1: np.ndarray = field(init=False)
    carried2: np.ndarray = field(init=False)
    nrec1: np.ndarray = field(init=False)
    nrec2: np.ndarray = field(init=False)
    an1: Optional[np.ndarray] = None
    an2: Optional[np.ndarray] = None
    anc_inconsistent: int = 0

    def __post_init__(self) -> None:
        n = self.n
        self.sum1 = np.zeros(n, np.int64)
        self.sum2 = np.zeros(n, np.int64)
        self.miss1 = np.zeros(n, bool)
        self.miss2 = np.zeros(n, bool)
        self.nonref_records = np.zeros(n, np.int16)
        self.unphased_hets = np.zeros(n, np.int16)
        # Bit 1: <50 bp gain, 2: <50 bp loss, 4: >=50 bp gain, 8: >=50 bp loss.
        self.carried1 = np.zeros(n, np.uint8)
        self.carried2 = np.zeros(n, np.uint8)
        self.nrec1 = np.zeros(n, np.int16)
        self.nrec2 = np.zeros(n, np.int16)

    def add(self, rec: Record) -> None:
        self.n_records += 1
        bp = rec.alt_bp
        cls = np.zeros(bp.size, np.uint8)
        small = np.abs(bp) < SV_MIN_BP
        cls[(bp > 0) & small] = 1
        cls[(bp < 0) & small] = 2
        cls[(bp > 0) & ~small] = 4
        cls[(bp < 0) & ~small] = 8
        for a, s, miss, carried, nrec in (
            (rec.a1, self.sum1, self.miss1, self.carried1, self.nrec1),
            (rec.a2, self.sum2, self.miss2, self.carried2, self.nrec2),
        ):
            ok = a >= 0
            idx = np.where(ok, a, 0)
            s += np.where(ok, bp[idx], 0)
            miss |= ~ok
            carried |= np.where(ok, cls[idx], 0).astype(np.uint8)
            nrec += (ok & (bp[idx] != 0)).astype(np.int16)
        nonref = ((rec.a1 > 0) & (bp[np.maximum(rec.a1, 0)] != 0)) | ((rec.a2 > 0) & (bp[np.maximum(rec.a2, 0)] != 0))
        self.nonref_records += nonref.astype(np.int16)
        both = (rec.a1 >= 0) & (rec.a2 >= 0)
        differ = both & (bp[np.maximum(rec.a1, 0)] != bp[np.maximum(rec.a2, 0)])
        self.unphased_hets += (rec.unphased_het & differ).astype(np.int16)
        if self.with_ancestry and rec.an1 is not None:
            if self.an1 is None:
                self.an1, self.an2 = rec.an1, rec.an2
            else:
                # A tandem array is orders of magnitude shorter than a FLARE
                # tract, so every record at a locus should carry the same local
                # ancestry. Count haplotypes where it does not instead of
                # assuming it: the first record's call is the one used.
                for cur, new in ((self.an1, rec.an1), (self.an2, rec.an2)):
                    both = (cur >= 0) & (new >= 0)
                    self.anc_inconsistent += int(np.count_nonzero(both & (cur != new)))


def open_text(path: str) -> TextIO:
    if path == "-":
        return sys.stdin
    return gzip.open(path, "rt") if path.endswith((".gz", ".bgz")) else open(path)


def parse_region(region: str) -> tuple[str, Optional[int], Optional[int]]:
    if ":" not in region:
        return region, None, None
    chrom, span = region.rsplit(":", 1)
    a, b = span.replace(",", "").split("-")
    return chrom, int(a), int(b)


def load_catalog(path: str, chrom: str, lo: Optional[int], hi: Optional[int]) -> list[Locus]:
    """Loci on chrom whose 1-based start falls in [lo, hi] (all when lo/hi are None)."""
    out: list[Locus] = []
    with open_text(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("#", "track", "browser")):
                continue
            f = line.rstrip("\n").split("\t")
            if f[0] != chrom:
                continue
            start, end = int(f[1]), int(f[2])
            if lo is not None and not (lo <= start + 1 <= hi):
                continue
            lid = f[3] if len(f) > 3 and f[3] else f"{chrom}:{start}-{end}"
            motifs = [m for m in (f[4] if len(f) > 4 else "").split(",") if m and m != "."]
            out.append(Locus(chrom, start, end, lid, motifs[0].upper() if motifs else "", len(motifs)))
    out.sort(key=lambda L: (L.start, L.end))
    return out


def _symbolic(alt: str) -> bool:
    return alt.startswith("<") or "[" in alt or "]" in alt


def allele_bp(ref: str, alts: Sequence[str], svtype: str, svlen: str, end: str, pos: int) -> Optional[np.ndarray]:
    """Signed bp change per allele index, or None when the record has no length change."""
    lens = [x for x in svlen.split(",")] if svlen not in ("", ".") else []
    st = svtype.split(",")[0].split(":")[0].upper() if svtype not in ("", ".") else ""
    out = [0]
    for i, alt in enumerate(alts):
        if alt in ("*", ".", "<*>", "<NON_REF>"):
            out.append(0)
            continue
        if not _symbolic(alt):
            out.append(len(alt) - len(ref))
            continue
        name = alt.strip("<>").split(":")[0].upper() if alt.startswith("<") else "BND"
        kind = name if name in ("DEL", "INS", "DUP", "CNV") else st
        raw = lens[i] if i < len(lens) else (lens[0] if len(lens) == 1 else "")
        size: Optional[int] = None
        if raw not in ("", "."):
            try:
                size = abs(int(float(raw)))
            except ValueError:
                size = None
        if size is None and kind == "DEL" and end not in ("", "."):
            size = max(int(end) - pos, 0)
        if size is None or kind not in ("DEL",) + GAIN_TYPES:
            out.append(0)
            continue
        out.append(-size if kind == "DEL" else size)
    arr = np.array(out, np.int64)
    return arr if np.any(arr[1:] != 0) else None


def record_span(pos: int, ref: str, alt_bp: np.ndarray, end: str) -> tuple[int, int]:
    start0 = pos - 1
    end0 = pos + 1
    if np.any(alt_bp < 0):
        e = start0 + len(ref)
        if end not in ("", "."):
            e = max(e, int(end))
        end0 = max(end0, e)
    return start0, end0


def _anc_codes(tokens: Sequence[str]) -> np.ndarray:
    """Ancestry codes as integers, -1 for '.' / empty. Any number of digits."""
    out = np.full(len(tokens), -1, np.int16)
    for i, t in enumerate(tokens):
        t = t.strip()
        if t and t != "." and (t.isdigit() or (t[0] == "-" and t[1:].isdigit())):
            out[i] = int(t)
    return out


def parse_gts(cells: str, n: int, with_ancestry: bool):
    """(a1, a2, unphased_het, an1, an2) for the per-sample columns of one query line.

    Ancestry codes are parsed as decimal integers of any width: a FLARE model
    with nanc >= 10 emits two-digit codes, and reading only the first byte would
    silently fold ancestry 12 into ancestry 1.
    """
    width = 8 if with_ancestry else 4
    if len(cells) == width * n - 1:
        b = np.frombuffer((cells + "\t").encode(), np.uint8).reshape(n, width)
        if np.all((b[:, 1] == 0x7C) | (b[:, 1] == 0x2F)) and np.all(b[:, 3] == 0x09):
            a1 = b[:, 0].astype(np.int16) - 48
            a2 = b[:, 2].astype(np.int16) - 48
            a1[(a1 < 0) | (a1 > 9)] = -1
            a2[(a2 < 0) | (a2 > 9)] = -1
            unph = (b[:, 1] == 0x2F) & (a1 != a2) & (a1 >= 0) & (a2 >= 0)
            an1 = an2 = None
            if with_ancestry:
                # Fixed-width path: every AN field is exactly one byte wide, so a
                # digit is its own code and anything else ('.') is missing.
                an1 = (b[:, 4].astype(np.int16) - 48)
                an2 = (b[:, 6].astype(np.int16) - 48)
                an1[(an1 < 0) | (an1 > 9)] = -1
                an2[(an2 < 0) | (an2 > 9)] = -1
            return a1, a2, unph, an1, an2
    parts = cells.split("\t")
    step = 3 if with_ancestry else 1
    if len(parts) != step * n:
        raise ValueError(f"expected {step * n} sample fields, got {len(parts)}")
    a1 = np.full(n, -1, np.int16)
    a2 = np.full(n, -1, np.int16)
    unph = np.zeros(n, bool)
    for i in range(n):
        gt = parts[step * i]
        sep = "|" if "|" in gt else "/"
        al = gt.split(sep)
        x = int(al[0]) if al[0].isdigit() else -1
        y = int(al[1]) if len(al) > 1 and al[1].isdigit() else -1
        a1[i], a2[i] = x, y
        unph[i] = sep == "/" and x != y and x >= 0 and y >= 0
    an1 = an2 = None
    if with_ancestry:
        an1 = _anc_codes(parts[1::3])
        an2 = _anc_codes(parts[2::3])
    return a1, a2, unph, an1, an2


class Source:
    """bcftools query stream for one VCF, reordered to the common sample order."""

    def __init__(self, vcf: str, samples_file: Optional[str], region_args: list[str], with_ancestry: bool,
                 with_ru: bool, order: Optional[list[str]] = None, prefilter: bool = True,
                 band: tuple[int, int] = (0, 1 << 62), filters: Sequence[str] = ()):
        self.vcf = vcf
        self.band = band
        fmt = "%CHROM\t%POS\t%ID\t%REF\t%ALT\t%INFO/SVTYPE\t%INFO/SVLEN\t%INFO/END"
        if with_ru:
            fmt += "\t%INFO/RU_TEST\t%INFO/PERIOD\t%INFO/MOTIF"
        fmt += "[\t%GT\t%AN1\t%AN2]\n" if with_ancestry else "[\t%GT]\n"
        # -u: a joint gVCF has no INFO/END, SVTYPE, or SVLEN. Without it, bcftools
        # exits instead of printing "." for those tags.
        cmd = ["bcftools", "query", "-u", "-H", "-f", fmt, *region_args]
        if samples_file:
            cmd += ["-S", samples_file, "--force-samples"]
        expr = []
        if prefilter:
            # gVCF reference blocks are ALT <NON_REF> / <*> and can carry a
            # megabase REF. Exclude them before strlen(), which would otherwise
            # keep every block (ALT~"<") and pull that REF into memory.
            expr.append('ALT!="<NON_REF>" && ALT!="<*>" && (strlen(REF)!=strlen(ALT) || ALT~"<")')
        if filters:
            expr.append("(" + " || ".join(f'FILTER="{f}"' for f in filters) + ")")
        if expr:
            cmd += ["-i", " && ".join(expr)]
        cmd.append(vcf)
        self.cmd = cmd
        self.n_fixed = 11 if with_ru else 8
        self.with_ancestry = with_ancestry
        self.proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
        assert self.proc.stdout is not None
        self.fh = self.proc.stdout
        header = self.fh.readline().rstrip("\n").lstrip("#").strip().split("\t")
        step = 3 if with_ancestry else 1
        names = [h.split("]", 1)[1].rsplit(":", 1)[0] for h in header[self.n_fixed::step]]
        self.samples = names
        self.perm: Optional[np.ndarray] = None
        if order is not None and names != order:
            idx = {s: i for i, s in enumerate(names)}
            self.perm = np.array([idx[s] for s in order], np.int64)

    def lines(self) -> Iterator[tuple[int, "Source", str]]:
        for line in self.fh:
            pos = int(line.split("\t", 2)[1])
            yield pos, self, line

    def close(self) -> None:
        self.fh.close()
        rc = self.proc.wait()
        if rc not in (0, -13):
            raise SystemExit(f"bcftools query failed ({rc}): {' '.join(self.cmd)}")


def vcf_samples(vcf: str) -> list[str]:
    out = subprocess.run(["bcftools", "query", "-l", vcf], check=True, capture_output=True, text=True).stdout
    return [s for s in out.split("\n") if s]


def contig_lines(vcf: str) -> list[str]:
    out = subprocess.run(["bcftools", "view", "-h", vcf], check=True, capture_output=True, text=True).stdout
    return [ln for ln in out.split("\n") if ln.startswith("##contig=")]


def fmt_units(v: float) -> str:
    """Human-readable units for the allele table."""
    return f"{v:.4g}" if v != int(v) else str(int(v))


def fmt_units_exact(v: float) -> str:
    """Units for --out-vcf INFO RU_DOSAGE, at round-trip precision.

    The association path consumes this value as the per-haplotype dosage, so it
    must not be quantised: at a period-3 locus a 7 bp change is 2.333333...
    units, and the allele table's 4 significant figures would introduce a
    ~1e-4 relative error into every test. SVLEN and PERIOD are also written, so
    a consumer can recompute the ratio exactly.
    """
    return f"{v:.17g}" if v != int(v) else str(int(v))


class Aggregator:
    def __init__(self, chrom: str, loci: list[Locus], n: int, *, flank: int, with_ancestry: bool,
                 missing_as_ref: bool, alleles_out: TextIO, vcf_out: Optional[TextIO]):
        self.chrom = chrom
        self.loci = loci
        self.starts = [L.start for L in loci]
        self.n = n
        self.flank = flank
        self.with_ancestry = with_ancestry
        self.missing_as_ref = missing_as_ref
        self.alleles_out = alleles_out
        self.vcf_out = vcf_out
        self.active: dict[int, LocusState] = {}
        self.pending: list[tuple[int, int, str, str]] = []
        self.counts: dict[str, int] = {
            k: 0 for k in (
                "records_read", "records_length_change", "records_assigned", "records_unassigned",
                "records_uncatalogued_ru", "loci_in_catalog", "loci_with_records", "loci_polymorphic",
                "loci_written", "hap_called", "hap_missing", "hap_unphased_ambiguous", "hap_multi_record",
                "hap_possible_duplicate", "alleles", "alleles_nonref", "alleles_out_of_band",
                "hap_ancestry_inconsistent", "hap_ancestry_missing", "anc_code_max",
            )
        }
        self.counts["anc_code_max"] = -1
        self.counts["loci_in_catalog"] = len(loci)

    def assign(self, start0: int, end0: int) -> Optional[int]:
        i = bisect.bisect_right(self.starts, end0 + self.flank) - 1
        best, best_key = None, None
        steps = 0
        while i >= 0 and steps < 8:
            L = self.loci[i]
            if L.end + self.flank >= start0:
                overlap = min(L.end, end0) - max(L.start, start0)
                dist = 0 if overlap > 0 else max(L.start - end0, start0 - L.end, 0)
                key = (overlap, -dist, -(L.end - L.start))
                if best_key is None or key > best_key:
                    best, best_key = i, key
            steps += 1
            i -= 1
        return best

    def flush_before(self, start0: int) -> None:
        done = [i for i, st in self.active.items() if st.locus.end + self.flank < start0]
        for i in done:
            self.finish(self.active.pop(i))
        self.emit_pending()

    def emit_pending(self, final: bool = False) -> None:
        if not self.pending:
            return
        floor = min((st.locus.start for st in self.active.values()), default=None)
        while self.pending and (final or floor is None or self.pending[0][0] <= floor):
            _start, _seq, allele_rows, vcf_line = heapq.heappop(self.pending)
            self.alleles_out.write(allele_rows)
            if vcf_line and self.vcf_out is not None:
                self.vcf_out.write(vcf_line)

    def add(self, rec: Record, idx: Optional[int]) -> None:
        if idx is None:
            if not rec.uncatalogued_period:
                self.counts["records_unassigned"] += 1
                return
            self.counts["records_uncatalogued_ru"] += 1
            motif = (rec.uncatalogued_motif or "").upper()
            L = Locus(self.chrom, rec.start0, rec.end0, rec.rid, motif, 1 if motif else 0, rec.uncatalogued_period)
            st = LocusState(L, self.n, self.with_ancestry)
            st.add(rec)
            self.finish(st)
            self.emit_pending()
            return
        self.counts["records_assigned"] += 1
        st = self.active.get(idx)
        if st is None:
            st = self.active[idx] = LocusState(self.loci[idx], self.n, self.with_ancestry)
        st.add(rec)

    def finish(self, st: LocusState) -> None:
        L = st.locus
        c = self.counts
        c["loci_with_records"] += 1
        c["hap_ancestry_inconsistent"] += st.anc_inconsistent
        if st.an1 is not None:
            for an in (st.an1, st.an2):
                c["hap_ancestry_missing"] += int(np.count_nonzero(an < 0))
                if an.size:
                    c["anc_code_max"] = max(c["anc_code_max"], int(an.max()))
        m1, m2 = st.miss1.copy(), st.miss2.copy()
        if self.missing_as_ref:
            m1[:] = False
            m2[:] = False
        ambiguous = st.unphased_hets >= 2
        m1 |= ambiguous
        m2 |= ambiguous
        c["hap_unphased_ambiguous"] += int(2 * ambiguous.sum())
        called1, called2 = ~m1, ~m2
        c["hap_multi_record"] += int(((st.nrec1 >= 2) & called1).sum() + ((st.nrec2 >= 2) & called2).sum())
        dup = lambda bits: ((bits & 1) & ((bits & 4) >> 2)) | (((bits & 2) >> 1) & ((bits & 8) >> 3))  # noqa: E731
        c["hap_possible_duplicate"] += int((dup(st.carried1) & called1).sum() + (dup(st.carried2) & called2).sum())
        vals = np.concatenate([st.sum1[called1], st.sum2[called2]])
        n_called = int(vals.size)
        n_missing = 2 * self.n - n_called
        c["hap_called"] += n_called
        c["hap_missing"] += n_missing
        uniq, cnt = np.unique(vals, return_counts=True) if n_called else (np.array([], np.int64), np.array([], np.int64))
        period = L.period
        head = f"{L.lid}\t{L.chrom}\t{L.start}\t{L.end}\t{period}\t{L.motif or '.'}\t{L.n_motifs}\t{st.n_records}\t{n_called}\t{n_missing}"
        rows = "".join(f"{head}\t{int(u)}\t{fmt_units(u / period)}\t{int(k)}\n" for u, k in zip(uniq, cnt))
        c["alleles"] += int(uniq.size)
        nonref = uniq[uniq != 0]
        c["alleles_nonref"] += int(nonref.size)
        vcf_line = ""
        if nonref.size:
            c["loci_polymorphic"] += 1
        if nonref.size and self.vcf_out is not None:
            index = {int(v): i + 1 for i, v in enumerate(nonref)}
            lut_keys = np.array(sorted(index), np.int64)
            def code(s: np.ndarray, called: np.ndarray) -> np.ndarray:
                out = np.zeros(s.size, np.int64)
                nz = s != 0
                out[nz] = np.searchsorted(lut_keys, s[nz]) + 1
                return np.where(called, out, -1)
            g1, g2 = code(st.sum1, called1), code(st.sum2, called2)
            alts = ",".join(f"<TR:{int(v):+d}bp>" for v in nonref)
            info = (
                f"END={L.end};RU_TEST;PERIOD={period};MOTIF={L.motif or '.'};LOCUS_RECORDS={st.n_records};"
                f"SVLEN={','.join(str(int(v)) for v in nonref)};"
                f"RU_DOSAGE={','.join(fmt_units_exact(v / period) for v in nonref)}"
            )
            gts = [f"{'.' if a < 0 else a}{'/' if amb else '|'}{'.' if b < 0 else b}" for a, b, amb in zip(g1, g2, ambiguous)]
            if self.with_ancestry and st.an1 is not None:
                fmt = "GT:AN1:AN2"
                cells = [
                    f"{g}:{'.' if x < 0 else int(x)}:{'.' if y < 0 else int(y)}"
                    for g, x, y in zip(gts, st.an1, st.an2)
                ]
            else:
                fmt = "GT"
                cells = gts
            vcf_line = "\t".join([L.chrom, str(L.start + 1), L.lid, "N", alts, ".", "PASS", info, fmt, *cells]) + "\n"
            c["loci_written"] += 1
        heapq.heappush(self.pending, (L.start, c["loci_with_records"], rows, vcf_line))

    def close(self) -> None:
        for i in sorted(self.active):
            self.finish(self.active[i])
        self.active.clear()
        self.emit_pending(final=True)


ALLELE_COLUMNS = (
    "locus_id", "chrom", "start", "end", "period", "motif", "n_motifs", "n_records",
    "n_hap", "n_hap_missing", "dosage_bp", "dosage_units", "n_hap_allele",
)


def vcf_header(contigs: list[str], samples: list[str], with_ancestry: bool) -> str:
    lines = ["##fileformat=VCFv4.2", '##FILTER=<ID=PASS,Description="All filters passed">', *contigs]
    lines += [
        '##ALT=<ID=TR,Description="Repeat locus allele; the suffix is the summed length change vs REF in bp">',
        '##INFO=<ID=END,Number=1,Type=Integer,Description="End of the catalog locus">',
        '##INFO=<ID=RU_TEST,Number=0,Type=Flag,Description="Repeat-unit dosage locus">',
        '##INFO=<ID=PERIOD,Number=1,Type=Integer,Description="Locus motif length in bp">',
        '##INFO=<ID=MOTIF,Number=1,Type=String,Description="Locus motif (first catalog motif)">',
        '##INFO=<ID=LOCUS_RECORDS,Number=1,Type=Integer,Description="Integrated-callset records summed into this locus">',
        '##INFO=<ID=SVLEN,Number=A,Type=Integer,Description="Summed length change vs REF in bp">',
        '##INFO=<ID=RU_DOSAGE,Number=A,Type=Float,Description="Signed repeat units vs REF (SVLEN / PERIOD)">',
        '##FORMAT=<ID=GT,Number=1,Type=String,Description="Locus allele per haplotype; unphased when assignment is ambiguous">',
    ]
    if with_ancestry:
        lines += [
            '##FORMAT=<ID=AN1,Number=1,Type=Integer,Description="FLARE ancestry, haplotype 1 (first record of the locus)">',
            '##FORMAT=<ID=AN2,Number=1,Type=Integer,Description="FLARE ancestry, haplotype 2 (first record of the locus)">',
        ]
    lines.append("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\t" + "\t".join(samples))
    return "\n".join(lines) + "\n"


def run(args: argparse.Namespace) -> dict:
    chrom, lo, hi = parse_region(args.region) if args.region else (args.chrom, None, None)
    if not chrom:
        raise SystemExit("--chrom or --region is required")
    loci = load_catalog(args.catalog_bed, chrom, lo, hi)
    print(f"[aggregate] {len(loci):,} catalog loci on {args.region or chrom}", file=sys.stderr, flush=True)

    keep: Optional[set[str]] = None
    if args.samples:
        keep = {s.strip() for s in Path(args.samples).read_text().split("\n") if s.strip()}
    prefixes = tuple(p for p in (args.exclude_prefixes or "").split(",") if p)
    split = args.split_bp
    inputs = (
        [(v, (0, 1 << 62)) for v in args.vcf or []]
        + [(v, (0, split)) for v in args.small_vcf or []]
        + [(v, (split, 1 << 62)) for v in args.sv_vcf or []]
    )
    if not inputs:
        raise SystemExit("give at least one --vcf, --small-vcf or --sv-vcf")
    filters = tuple(f for f in (args.apply_filters or "").split(",") if f)
    per_vcf = [vcf_samples(v) for v, _b in inputs]
    common = [s for s in per_vcf[0] if all(s in set(o) for o in per_vcf[1:])]
    if keep is not None:
        common = [s for s in common if s in keep]
    if prefixes:
        common = [s for s in common if not s.startswith(prefixes)]
    if not common:
        raise SystemExit("no samples left after --samples / --exclude-prefixes")
    n = len(common)

    tmp = tempfile.TemporaryDirectory()
    samples_file = os.path.join(tmp.name, "samples.txt")
    Path(samples_file).write_text("\n".join(common) + "\n")
    # One index jump to the window (each record once, by POS); unless uncatalogued RU_TEST
    # records are kept, stream only records overlapping a catalog locus. -T streams, whereas
    # -R with one region per locus seeks per locus and is several times slower.
    span = chrom
    if lo is not None:
        first = min((L.start for L in loci), default=lo)
        last = max((L.end for L in loci), default=hi)
        # Window-start pad reaches deletions that begin before the first locus.
        span = f"{chrom}:{max(first - args.window_pad, 1)}-{last + args.flank}"
    region_args = ["-r", span, "--regions-overlap", "pos"]
    if loci and not args.keep_uncatalogued_ru:
        bed = os.path.join(tmp.name, "loci.bed")
        with open(bed, "w") as fh:
            for L in loci:
                fh.write(f"{chrom}\t{max(L.start - args.flank - 1, 0)}\t{L.end + args.flank}\n")
        region_args += ["-T", bed, "--targets-overlap", "variant"]

    alleles_out = gzip.open(args.out_alleles, "wt") if args.out_alleles.endswith(".gz") else open(args.out_alleles, "w")
    alleles_out.write("\t".join(ALLELE_COLUMNS) + "\n")
    vcf_out = None
    if args.out_vcf:
        vcf_out = gzip.open(args.out_vcf, "wt") if args.out_vcf.endswith(".gz") and not args.out_vcf.endswith(".bgz") else open(args.out_vcf, "w")
        vcf_out.write(vcf_header(contig_lines(inputs[0][0]), common, args.with_ancestry))

    agg = Aggregator(chrom, loci, n, flank=args.flank, with_ancestry=args.with_ancestry,
                     missing_as_ref=args.missing_as_ref, alleles_out=alleles_out, vcf_out=vcf_out)
    sources = [
        Source(v, samples_file, region_args, args.with_ancestry, args.keep_uncatalogued_ru, order=common,
               band=band, filters=filters)
        for v, band in inputs
    ]
    seen: set[tuple[int, str, str, str, tuple[int, int]]] = set()
    last_pos = -1
    for pos, src, line in heapq.merge(*(s.lines() for s in sources), key=lambda t: t[0]):
        agg.counts["records_read"] += 1
        src_fixed = sources[0].n_fixed
        f = line.rstrip("\n").split("\t", src_fixed)
        rid, ref, alt_s, svtype, svlen, end = f[2], f[3], f[4], f[5], f[6], f[7]
        if pos != last_pos:
            seen = set()
            last_pos = pos
        # A record in both a small and an SV source is kept from each; the size band keeps one allele.
        key = (pos, ref, alt_s, rid, src.band)
        if key in seen:
            continue
        seen.add(key)
        alts = alt_s.split(",")
        bp = allele_bp(ref, alts, svtype, svlen, end, pos)
        if bp is None or np.abs(bp).max() > args.max_bp:
            continue
        lo_bp, hi_bp = src.band
        size = np.abs(bp)
        out_band = (bp != 0) & ((size < lo_bp) | (size >= hi_bp))
        if out_band.any():
            agg.counts["alleles_out_of_band"] += int(out_band.sum())
            bp = np.where(out_band, 0, bp)
            if not bp.any():
                continue
        agg.counts["records_length_change"] += 1
        start0, end0 = record_span(pos, ref, bp, end)
        agg.flush_before(start0)
        cells = f[src_fixed] if len(f) > src_fixed else ""
        a1, a2, unph, an1, an2 = parse_gts(cells, n, args.with_ancestry)
        if args.num_ancs and an1 is not None:
            for an in (an1, an2):
                if an.size and an.max() >= args.num_ancs:
                    raise SystemExit(
                        f"{chrom}:{pos}: ancestry code {int(an.max())} but --num-ancs={args.num_ancs}"
                    )
        if args.ignore_phase:
            unph = (a1 != a2) & (a1 >= 0) & (a2 >= 0)
        if a1.size and (a1.max() >= bp.size or a2.max() >= bp.size):
            raise SystemExit(f"{chrom}:{pos}: allele index beyond ALT list")
        rec = Record(pos, start0, end0, rid if rid != "." else f"{chrom}:{pos}", bp, a1, a2, unph, an1, an2)
        idx = agg.assign(start0, end0) if loci else None
        if idx is None and args.keep_uncatalogued_ru:
            ru_test, period = f[8], f[9]
            if ru_test not in (".", "0", "") and period not in (".", ""):
                rec.uncatalogued_period = int(period)
                rec.uncatalogued_motif = f[10] if f[10] not in (".", "") else None
        agg.add(rec, idx)
    for s in sources:
        s.close()
    agg.close()
    alleles_out.close()
    if vcf_out is not None:
        vcf_out.close()
    tmp.cleanup()

    c = agg.counts
    summary = {
        "region": args.region or chrom,
        "vcfs": {"all": args.vcf or [], "small": args.small_vcf or [], "sv": args.sv_vcf or []},
        "split_bp": split,
        "apply_filters": list(filters),
        "catalog_bed": args.catalog_bed,
        "n_samples": n,
        "flank": args.flank,
        "with_ancestry": args.with_ancestry,
        "missing_as_ref": args.missing_as_ref,
        "ignore_phase": args.ignore_phase,
        "counts": c,
        "hap_missing_frac": c["hap_missing"] / max(c["hap_called"] + c["hap_missing"], 1),
        "hap_possible_duplicate_frac": c["hap_possible_duplicate"] / max(c["hap_called"], 1),
        # Fraction of haplotypes where a later record at the locus disagreed with
        # the first record's local-ancestry call. Should be ~0: a tandem array is
        # far shorter than a FLARE tract. A non-trivial value means the
        # "ancestry is constant across the locus" assumption is not holding.
        "hap_ancestry_inconsistent_frac": (
            c["hap_ancestry_inconsistent"] / max(c["hap_called"] + c["hap_missing"], 1)
        ),
        "anc_code_max": c["anc_code_max"],
    }
    return summary


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", action="append", help="Indexed VCF/BCF; every length change counts (repeatable)")
    p.add_argument("--small-vcf", action="append",
                   help="Small-variant VCF; only alleles with |change| < --split-bp count (repeatable)")
    p.add_argument("--sv-vcf", action="append",
                   help="SV VCF/BCF; only alleles with |change| >= --split-bp count (repeatable)")
    p.add_argument("--split-bp", type=int, default=50,
                   help="Size boundary between --small-vcf and --sv-vcf, so an event is counted from one callset")
    p.add_argument("--apply-filters", default="", help="Keep records whose FILTER is one of these, e.g. PASS,.")
    p.add_argument("--catalog-bed", required=True, help="chrom start end id motifs (BED, .gz ok)")
    p.add_argument("--chrom", default=None)
    p.add_argument("--region", default=None, help="chrom:start-end; loci whose start falls here")
    p.add_argument("--samples", default=None, help="Keep-list of sample ids")
    p.add_argument("--exclude-prefixes", default="", help="Drop samples with these id prefixes, e.g. HG,NA")
    p.add_argument("--flank", type=int, default=10, help="bp around a locus that still assigns a record to it")
    p.add_argument("--max-bp", type=int, default=100_000, help="Ignore records with a larger length change")
    p.add_argument("--window-pad", type=int, default=1000,
                   help="With --region, also read records starting this far before the first locus")
    p.add_argument("--with-ancestry", action="store_true", help="Carry FORMAT AN1/AN2 into --out-vcf")
    p.add_argument("--num-ancs", type=int, default=0,
                   help="With --with-ancestry, fail if any AN1/AN2 code is >= this (0 = no check)")
    p.add_argument("--keep-uncatalogued-ru", action="store_true",
                   help="Write RU_TEST records outside every catalog locus as their own loci")
    p.add_argument("--missing-as-ref", action="store_true", help="Treat missing calls as REF instead of missing")
    p.add_argument("--ignore-phase", action="store_true",
                   help="Treat phased hets as unphased, to score a phased callset like an unphased one")
    p.add_argument("--out-alleles", required=True)
    p.add_argument("--out-vcf", default=None)
    p.add_argument("--out-summary", required=True)
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    summary = run(args)
    Path(args.out_summary).write_text(json.dumps(summary, indent=2) + "\n")
    print(json.dumps(summary["counts"]), file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
