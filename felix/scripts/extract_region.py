#!/usr/bin/env python3
"""Pull one genomic window out of a FelixGenome run's merged results, for each arm, without downloading the files.

The merged per-phenotype results are 0.5-27 GB (merged.felix.tsv, or merged.felix.tsv.gz for runs since 2026-10-15), so this
streams ``gsutil cat`` and keeps only the rows whose CHR and POS fall in the window. Written for the Figure 3 locus card (regional plot, per-ancestry effects), but it takes any window.

    python3 extract_region.py --root gs://BUCKET/submissions/SUB/FelixGenome/WF \
        --chrom 9 --start 8480000 --end 9480000 --out-prefix fig3_locus_BI_174

``--root`` is the workflow directory; take it from the run's run_summary_tar URI by cutting everything from /call-Bundle.
Writes ``<out-prefix>.main.tsv``, ``.ru.tsv`` and ``.ru_biallelic.tsv`` (header plus the rows in the window, every column).
An arm whose merged file is absent is reported and skipped.
"""

from __future__ import annotations

import argparse
import gzip
import io
import re
import subprocess
import sys
from pathlib import Path
from typing import Iterable, Iterator

ARM_CALL = {"main": "Concat", "ru": "RuConcat", "ru_biallelic": "BaConcat"}
MERGED = "merged.felix.tsv"
# Since 712c4a2 (2026-10-15) the concat tasks write merged.felix.tsv.gz; runs before that wrote the plain file. Both are read.
MERGED_NAMES = (MERGED + ".gz", MERGED)


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


def list_call(root: str, call: str) -> list[tuple[int, str]]:
    """(size, uri) of every object under <root>/call-<call>/ (gsutil ls -l -r); empty when nothing is listed."""
    out = subprocess.run(["gsutil", "ls", "-l", "-r", f"{root.rstrip('/')}/call-{call}/"], capture_output=True, text=True)
    objs = []
    for line in out.stdout.splitlines():
        parts = line.split()
        if len(parts) >= 3 and parts[0].isdigit() and parts[-1].startswith("gs://"):
            objs.append((int(parts[0]), parts[-1]))
    return objs


def find_merged(objs: list[tuple[int, str]]) -> str | None:
    """Largest merged table (plain or .gz) among the listed objects: the real file, not a stub or a smaller attempt."""
    best = None
    for size, uri in objs:
        if any(uri.endswith("/" + n) for n in MERGED_NAMES) and (best is None or size > best[0]):
            best = (size, uri)
    return best[1] if best else None


def extract(uri: str, out: Path, chrom: str, start: int, end: int) -> int:
    cat = subprocess.Popen(["gsutil", "cat", uri], stdout=subprocess.PIPE, bufsize=1 << 20)
    stream = gzip.GzipFile(fileobj=cat.stdout) if uri.endswith(".gz") else cat.stdout
    n = -1  # the header is not a row
    with open(out, "w") as fh:
        for line in filter_region(io.TextIOWrapper(stream, encoding="utf-8"), chrom, start, end):
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
        objs = list_call(a.root, ARM_CALL[arm])
        uri = find_merged(objs)
        if not uri:
            print(f"{arm}: no {' or '.join(MERGED_NAMES)} under call-{ARM_CALL[arm]}/ in {a.root}", file=sys.stderr)
            print("   gsutil listed " + (f"{len(objs)} objects there, e.g. " + ", ".join(u.rsplit('/', 1)[-1] for _, u in objs[:6])
                                         if objs else "nothing there: check the root and your access"), file=sys.stderr)
            bad += 1
            continue
        out = Path(f"{a.out_prefix}.{arm}.tsv")
        n = extract(uri, out, a.chrom, a.start, a.end)
        print(f"{arm}: {n:,} rows in chr{a.chrom}:{a.start:,}-{a.end:,} -> {out}  (from {uri})")
    return 1 if bad else 0


if __name__ == "__main__":
    raise SystemExit(main())
