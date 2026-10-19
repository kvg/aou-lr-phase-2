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


def split_windows(lines: Iterable[str], windows: list[tuple[str, int, int]], chr_col="CHR", pos_col="POS"):
    """Yield (window index, line) for every data line falling in a window, header first as index -1.

    One pass over the stream serves every window, which is the point: a run's merged table is 0.5-27 GB, so reading it
    once for all of that phenotype's loci is the difference between minutes and an hour. A line inside two overlapping
    windows is yielded once per window.
    """
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
    want = [(_chrom_key(str(c)), s, e) for c, s, e in windows]
    yield -1, header
    for line in it:
        f = line.split("\t", upto)
        if len(f) <= upto - 1:
            continue
        key = _chrom_key(f[ic])
        if not any(key == c for c, _, _ in want):
            continue
        try:
            pos = int(float(f[ip]))
        except ValueError:
            continue
        for i, (c, s, e) in enumerate(want):
            if key == c and s <= pos <= e:
                yield i, line


def _open_stream(uri: str):
    """(process, text stream) for a gs:// object or a local file; the caller closes the process."""
    if not str(uri).startswith("gs://"):
        fh = gzip.open(uri, "rt") if str(uri).endswith(".gz") else open(uri)
        return None, fh
    cat = subprocess.Popen(["gsutil", "cat", uri], stdout=subprocess.PIPE, bufsize=1 << 20)
    stream = gzip.GzipFile(fileobj=cat.stdout) if uri.endswith(".gz") else cat.stdout
    return cat, io.TextIOWrapper(stream, encoding="utf-8")


def extract_many(uri: str, targets: list[tuple[str, int, int, Path]]) -> list[int]:
    """Cut several windows out of one table in a single pass; returns the row count written per window."""
    windows = [(c, s, e) for c, s, e, _ in targets]
    cat, stream = _open_stream(uri)
    handles = [open(p, "w") for *_, p in targets]
    counts = [0] * len(targets)
    try:
        for idx, line in split_windows(stream, windows):
            if idx < 0:
                for fh in handles:
                    fh.write(line)
            else:
                handles[idx].write(line)
                counts[idx] += 1
    finally:
        for fh in handles:
            fh.close()
        if cat is not None:
            cat.stdout.close()
            if cat.wait() != 0:
                raise RuntimeError(f"gsutil cat failed for {uri}")
        else:
            stream.close()
    return counts


def extract(uri: str, out: Path, chrom: str, start: int, end: int) -> int:
    """One window; kept as the single-window entry point used by the Figure 3 card cell."""
    return extract_many(uri, [(chrom, start, end, out)])[0]


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--root", required=True, help="gs://.../submissions/<sub>/FelixGenome/<workflow>")
    ap.add_argument("--chrom")
    ap.add_argument("--start", type=int)
    ap.add_argument("--end", type=int)
    ap.add_argument("--out-prefix")
    ap.add_argument("--windows", help="TSV with chrom, start, end, out_prefix: several windows cut in ONE pass "
                                      "per arm (use instead of --chrom/--start/--end/--out-prefix)")
    ap.add_argument("--arms", default=",".join(ARM_CALL), help="comma-separated subset of main,ru,ru_biallelic")
    a = ap.parse_args()

    if a.windows:
        targets = read_windows(Path(a.windows))
        if not targets:
            ap.error(f"{a.windows} lists no windows")
    else:
        if not all([a.chrom, a.start is not None, a.end is not None, a.out_prefix]):
            ap.error("give --windows, or all of --chrom --start --end --out-prefix")
        targets = [(a.chrom, a.start, a.end, a.out_prefix)]

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
        outs = [(c, s, e, Path(f"{pre}.{arm}.tsv")) for c, s, e, pre in targets]
        counts = extract_many(uri, outs)
        for (c, s, e, out), n in zip(outs, counts):
            print(f"{arm}: {n:,} rows in chr{c}:{s:,}-{e:,} -> {out}  (from {uri})")
    return 1 if bad else 0


def read_windows(path: Path) -> list[tuple[str, int, int, str]]:
    """(chrom, start, end, out_prefix) per line of a windows TSV; '#' comments and a header row are skipped."""
    out = []
    for line in path.read_text().splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        f = line.split("\t")
        if len(f) < 4:
            raise ValueError(f"{path}: need 4 tab-separated fields (chrom, start, end, out_prefix), got {line!r}")
        if not f[1].strip().isdigit():  # the header
            continue
        out.append((f[0].strip(), int(f[1]), int(f[2]), f[3].strip()))
    return out


if __name__ == "__main__":
    raise SystemExit(main())
