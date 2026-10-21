#!/usr/bin/env python3
"""Rebuild TRGT's repeat-definition BED from a TRGT output VCF.

``trgt --repeats`` needs four tab-separated fields::

    chrom  start  end  ID=<trid>;MOTIFS=<motifs>;STRUC=<struc>

The catalog in ``refs/repeat_loci/`` is five fields (``chrom start end TRID
motifs``) because ``trgt_plvi.py catalog`` wrote it for
``aggregate_repeat_loci.py``, dropping ``STRUC`` on the way. ``STRUC`` cannot be
reconstructed from a motif list at a compound locus, so that file cannot be
converted — but a TRGT *output* VCF carries all three fields in INFO, because
TRGT echoes the definition it was given::

    TRID=1-66431-66446-TATAT;END=66446;MOTIFS=TATAT;STRUC=<TR22>
    TRID=1-66160-66175-TA,1-66192-66201-TA,...;END=66415;MOTIFS=ATATA,TATAT,TA;STRUC=<VC2>

So one genotyped sample round-trips the definitions it was genotyped against,
which is also the set that makes a re-genotyping comparison locus-for-locus.

    trgt_defs_from_vcf.py --vcf <sample>_trgt.TR_Explorer_1.01.vcf.gz \\
        --out trexplorer_v1.0.1.trgt_definitions.bed

A ``gs://`` input is copied down with gsutil first. htslib in the Terra
notebook image cannot open ``gs://`` at all -- it reports "Protocol not
supported" even with GCS_OAUTH_TOKEN and GCS_REQUESTER_PAYS_PROJECT exported,
so that is a build limitation, not a credentials one. A TRGT VCF is ~180 MB.

Coordinates: VCF ``POS`` is the 1-based locus start and BED is half-open, so the
start is ``POS - 1`` and the end is INFO ``END``. This matches the derived
catalog, whose rows read ``chr20 10000929 10000950 20-10000930-10000950-A``.

**One sample's VCF holds the loci TRGT emitted for that sample.** Any locus it
dropped is absent here. Check the row count against the catalog ceiling
(``gsutil cat <catalog> | gunzip | wc -l``) before treating the output as the
whole catalog; ``--expect-loci`` turns that check into a hard one.
"""

from __future__ import annotations

import argparse
import gzip
import os
import re
import shutil
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Iterator, Optional

QUERY = "%CHROM\t%POS\t%INFO/END\t%INFO/TRID\t%INFO/MOTIFS\t%INFO/STRUC\n"
MISSING = (".", "")


def query_lines(vcf: str) -> Iterator[str]:
    """bcftools query stream. stderr is surfaced: an unreadable URI must not look empty."""
    proc = subprocess.Popen(["bcftools", "query", "-f", QUERY, vcf],
                            stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True,
                            bufsize=1 << 20)
    assert proc.stdout is not None and proc.stderr is not None
    try:
        yield from proc.stdout
    finally:
        proc.stdout.close()
        err = proc.stderr.read()
        proc.stderr.close()
        if proc.wait() != 0:
            raise SystemExit(f"bcftools query failed on {vcf}:\n{err.strip()}")


def localize(uri: str, into: str) -> str:
    """Copy a gs:// object down and return the local path.

    htslib in the Terra notebook image cannot open gs:// at all -- it reports
    `Protocol not supported` even with GCS_OAUTH_TOKEN and
    GCS_REQUESTER_PAYS_PROJECT exported, so this is a build limitation and not a
    credentials one. gsutil works, which is why the rest of this repo streams
    with `gsutil cat` rather than through htslib. A TRGT VCF is ~180 MB, so
    copying it is cheaper than arguing with the protocol.
    """
    dest = os.path.join(into, os.path.basename(uri))
    cmd = ["gsutil", "-q", "cp", uri, dest]
    proc = subprocess.run(cmd, capture_output=True, text=True)
    if proc.returncode != 0 or not os.path.exists(dest):
        raise SystemExit(f"gsutil cp failed for {uri}:\n{proc.stderr.strip()}")
    # bcftools wants the index beside the file for random access; this script
    # reads sequentially, so a missing .tbi is not an error.
    for suffix in (".tbi", ".csi"):
        idx = subprocess.run(["gsutil", "-q", "cp", uri + suffix, dest + suffix],
                             capture_output=True, text=True)
        if idx.returncode == 0:
            break
    return dest


TRID_COORDS = re.compile(r"^[^-]+-(\d+)-(\d+)-")


def definition(chrom: str, pos: str, end: str, trid: str, motifs: str, struc: str) -> str:
    """One TRGT repeat-definition row.

    The start is POS unchanged, NOT POS - 1. TRGT writes the BED start into
    the VCF POS verbatim, so a round trip is only the identity when the two
    are equal. This was measured, not read off a spec: feeding definitions
    built with POS - 1 shifted every locus 1 bp left, which showed up as a
    POS delta of +1 at all 14,618 loci in a chr20 window against the
    production call, with END unchanged and AL differing at 14,613 of them.
    The repeat framing moves with the anchor, so the error is not confined
    to the reference span -- and under REF-relative coding a uniform shift
    adds c * ANC{k} to every per-ancestry dosage, the defect REPEAT_DOSAGE.md
    section 2.1 exists to prevent.
    """
    return f"{chrom}\t{int(pos)}\t{int(end)}\tID={trid};MOTIFS={motifs};STRUC={struc}"


def trid_coords(trid: str) -> Optional[tuple[int, int]]:
    """Span encoded in a TRID, as (start of the first part, end of the last).

    A compound locus carries one ``chrom-start-end-motif`` part per component,
    comma-joined, so the locus span runs from the FIRST part's start to the
    LAST part's end -- not the first part's end.

    This is a weak check by design. The span is not an exact invariant of the
    TRExplorer catalog: some loci begin tens of bp before their first component
    and some single-component loci end after their TRID says. It is here only
    to catch a SYSTEMATIC shift, which drives agreement to ~0 while those
    idiosyncrasies leave it high -- hence the rate threshold in main() rather
    than a per-row assertion.
    """
    parts = trid.split(",")
    first, last = TRID_COORDS.match(parts[0]), TRID_COORDS.match(parts[-1])
    if first is None or last is None:
        return None
    return int(first.group(1)), int(last.group(2))


def run(args: argparse.Namespace) -> dict:
    counts = {"records": 0, "written": 0, "skipped_missing_field": 0,
              "trid_checked": 0, "trid_start_agree": 0, "trid_span_agree": 0}
    skipped: list[str] = []
    mismatches: list[str] = []
    scratch = tempfile.mkdtemp(prefix="trgt_defs_")
    try:
        vcf = localize(args.vcf, scratch) if args.vcf.startswith("gs://") else args.vcf
    except BaseException:
        shutil.rmtree(scratch, ignore_errors=True)
        raise
    # Write to a sibling temp file and rename only on success. The first
    # version opened args.out directly, so a bcftools failure left a valid
    # but empty .gz behind -- which was then uploaded over the real path.
    tmp = args.out + ".tmp"
    opener = gzip.open if args.out.endswith(".gz") else open
    try:
        with opener(tmp, "wt") as fh:
            for line in query_lines(vcf):
                fields = line.rstrip("\n").split("\t")
                if len(fields) < 6:
                    continue
                counts["records"] += 1
                chrom, pos, end, trid, motifs, struc = fields[:6]
                if any(v in MISSING for v in (end, trid, motifs, struc)):
                    counts["skipped_missing_field"] += 1
                    if len(skipped) < 5:
                        skipped.append(f"{chrom}:{pos} TRID={trid} MOTIFS={motifs} STRUC={struc}")
                    continue
                row = definition(chrom, pos, end, trid, motifs, struc)
                coords = trid_coords(trid)
                if coords is not None:
                    counts["trid_checked"] += 1
                    fields = row.split("\t")
                    got = (int(fields[1]), int(fields[2]))
                    if got[0] == coords[0]:
                        counts["trid_start_agree"] += 1
                    elif len(mismatches) < 5:
                        mismatches.append(
                            f"{trid.split(',')[0]}...: wrote start {got[0]}, "
                            f"TRID says {coords[0]}")
                    if got == coords:
                        counts["trid_span_agree"] += 1
                fh.write(row + "\n")
                counts["written"] += 1
    except BaseException:
        # Leave nothing behind. A partial file at args.out is worse than none:
        # the first version's empty .gz got uploaded over the canonical path.
        if os.path.exists(tmp):
            os.unlink(tmp)
        raise
    finally:
        shutil.rmtree(scratch, ignore_errors=True)

    checked = counts["trid_checked"]
    if checked:
        rate = counts["trid_start_agree"] / checked
        counts["trid_start_agree_frac"] = round(rate, 6)
        counts["trid_span_agree_frac"] = round(counts["trid_span_agree"] / checked, 6)
        if rate < args.min_trid_agreement:
            os.unlink(tmp)
            raise SystemExit(
                f"error: the start agrees with the locus id on only {rate:.1%} of "
                f"{checked:,} rows, below the --min-trid-agreement floor of "
                f"{args.min_trid_agreement:.0%}. A systematic coordinate shift drives "
                f"this to ~0; catalog idiosyncrasies leave it high. If the rate "
                f"above is high but under the floor, that is catalog structure -- "
                f"re-run with --min-trid-agreement 0 to write anyway. Examples:\n  "
                + "\n  ".join(mismatches))
    if counts["written"] == 0:
        os.unlink(tmp)
        raise SystemExit(
            f"error: no definition written from {args.vcf}. If every record lacked STRUC, this "
            "is not a TRGT output VCF, or TRGT was run without a STRUC-bearing catalog.")
    if counts["skipped_missing_field"]:
        print(f"warning: skipped {counts['skipped_missing_field']:,} record(s) missing END, TRID, "
              f"MOTIFS or STRUC; e.g. {'; '.join(skipped)}", file=sys.stderr)
    if args.expect_loci is not None and counts["written"] != args.expect_loci:
        os.unlink(tmp)
        raise SystemExit(
            f"error: wrote {counts['written']:,} definitions but --expect-loci is "
            f"{args.expect_loci:,}. One sample's VCF holds only the loci TRGT emitted for it.")

    os.replace(tmp, args.out)
    return counts


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--vcf", required=True, help="A TRGT output VCF (local path or gs:// URI)")
    p.add_argument("--out", required=True, help="Output BED; .gz is compressed")
    p.add_argument(
        "--min-trid-agreement", type=float, default=0.90,
        help="Floor on the fraction of rows whose start matches the coordinate "
             "encoded in their own TRID. Catches a systematic coordinate shift "
             "(which drives it to ~0) without tripping on the few per cent of "
             "TRExplorer loci whose TRID span is idiosyncratic. Default 0.90.")
    p.add_argument("--expect-loci", type=int, default=None,
                   help="Fail unless exactly this many definitions are written, e.g. the "
                        "catalog's line count")
    return p


def main(argv: Optional[list[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    counts = run(args)
    for k, v in counts.items():
        print(f"{k:>24}  {v:,}")
    print(f"wrote {args.out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
