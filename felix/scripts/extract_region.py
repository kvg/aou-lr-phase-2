#!/usr/bin/env python3
"""Pull one genomic window out of a FelixGenome run's merged results, for each arm, without downloading the files.

The merged per-phenotype results are 0.5-27 GB, so this streams ``gsutil cat`` and keeps only the rows whose CHR and POS
fall in the window. Written for the Figure 3 locus card (regional plot, per-ancestry effects), but it takes any window.

    python3 extract_region.py --root gs://BUCKET/submissions/SUB/FelixGenome/WF \
        --chrom 9 --start 8480000 --end 9480000 --out-prefix fig3_locus_BI_174

``--root`` is the workflow directory; take it from the run's run_summary_tar URI by cutting everything from /call-Bundle.
Writes ``<out-prefix>.main.tsv``, ``.ru.tsv`` and ``.ru_biallelic.tsv`` (header plus the rows in the window, every column).
An arm whose merged file is absent is reported and skipped.
"""

from __future__ import annotations

import argparse
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Iterator

ARM_CALL = {"main": "Concat", "ru": "RuConcat", "ru_biallelic": "BaConcat"}
MERGED = "merged.felix.tsv"


def _chrom_key(v: str) -> str:
    return re.sub(r"^chr", "", v.strip(), flags=re.IGNORECASE)


def filter_region(lines: Iterable[str], chrom: str, start: int, end: int, chr_col="CHR", pos_col="POS") -> Iterator[str]:
    """Yield the header, then every data line with CHR == chrom and start <= POS <= end. Lines keep their newline."""
    it = iter(lines)
    header = next(it, None)
    if header is None:
        return
    cols = header.rstrip("\n").split("\t")
    for need in (chr_col, pos_col):
        if need not in cols:
            raise KeyError(f"no {need} column in the header: {cols[:12]}...")
    ic, ip = cols.index(chr_col), cols.index(pos_col)
    upto = max(ic, ip) + 1
    want = _chrom_key(str(chrom))
    yield header
    for line in it:
        f = line.split("\t", upto)  # split only as far as the two columns we need
        if len(f) <= upto - 1:
            continue
        if _chrom_key(f[ic]) != want:
            continue
        try:
            pos = int(float(f[ip]))
        except ValueError:
            continue
        if start <= pos <= end:
            yield line


def find_merged(root: str, call: str) -> str | None:
    """Largest merged.felix.tsv under <root>/call-<call>/ (the real file, not a stub or a smaller attempt)."""
    out = subprocess.run(["gsutil", "ls", "-l", "-r", f"{root.rstrip('/')}/call-{call}/"], capture_output=True, text=True)
    best = None
    for line in out.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[-1].endswith("/" + MERGED) and parts[0].isdigit():
            if best is None or int(parts[0]) > best[0]:
                best = (int(parts[0]), parts[-1])
    return best[1] if best else None


def extract(uri: str, out: Path, chrom: str, start: int, end: int) -> int:
    cat = subprocess.Popen(["gsutil", "cat", uri], stdout=subprocess.PIPE, text=True, bufsize=1 << 20)
    n = -1  # the header is not a row
    with open(out, "w") as fh:
        for line in filter_region(cat.stdout, chrom, start, end):
            fh.write(line)
            n += 1
    cat.stdout.close()
    if cat.wait() != 0:
        raise RuntimeError(f"gsutil cat failed for {uri}")
    return max(n, 0)


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="gs://.../submissions/<sub>/FelixGenome/<workflow>")
    ap.add_argument("--chrom", required=True)
    ap.add_argument("--start", type=int, required=True)
    ap.add_argument("--end", type=int, required=True)
    ap.add_argument("--out-prefix", required=True)
    ap.add_argument("--arms", default=",".join(ARM_CALL), help="comma-separated subset of main,ru,ru_biallelic")
    a = ap.parse_args()
    bad = 0
    for arm in a.arms.split(","):
        uri = find_merged(a.root, ARM_CALL[arm])
        if not uri:
            print(f"{arm}: no {MERGED} under call-{ARM_CALL[arm]}/ in {a.root}", file=sys.stderr)
            bad += 1
            continue
        out = Path(f"{a.out_prefix}.{arm}.tsv")
        n = extract(uri, out, a.chrom, a.start, a.end)
        print(f"{arm}: {n:,} rows in chr{a.chrom}:{a.start:,}-{a.end:,} -> {out}  (from {uri})")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
