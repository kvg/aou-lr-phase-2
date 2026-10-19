"""Tests for extract_region.py: the window filter, the choice of merged file, and streaming plain and gzip tables (gsutil is faked)."""

import gzip
import io

import pytest

import extract_region as er

HEADER = "MarkerID\tCHR\tPOS\tP_cct_admixed_c\n"
ROWS = [
    "a\tchr9\t8479999\t0.5\n",   # one base before the window
    "b\tchr9\t8480000\t0.4\n",   # first base of the window
    "c\tchr9\t8991058\t1e-9\n",
    "d\tchr9\t9480000\t0.3\n",   # last base of the window
    "e\tchr9\t9480001\t0.2\n",   # one base after
    "f\tchr19\t8991058\t0.1\n",  # right position, wrong chromosome (19 must not match 9)
    "g\tchr10\t8991058\t0.1\n",
]


def kept(lines, chrom="9", start=8_480_000, end=9_480_000):
    return [x.split("\t")[0] for x in list(er.filter_region(lines, chrom, start, end))[1:]]


def test_window_is_inclusive_and_chromosome_exact():
    assert kept([HEADER] + ROWS) == ["b", "c", "d"]


def test_chr_prefix_does_not_matter():
    bare = [r.replace("chr", "") for r in ROWS]
    assert kept([HEADER] + bare) == ["b", "c", "d"]
    assert kept([HEADER] + ROWS, chrom="chr9") == ["b", "c", "d"]


def test_header_is_always_first_and_unchanged():
    assert list(er.filter_region([HEADER] + ROWS, "9", 0, 10))[0] == HEADER
    assert list(er.filter_region([HEADER], "9", 0, 10)) == [HEADER]


def test_columns_after_the_two_we_need_are_kept_intact():
    wide = "MarkerID\tCHR\tPOS\t" + "\t".join(f"c{i}" for i in range(150)) + "\n"
    row = "x\t9\t100\t" + "\t".join(str(i) for i in range(150)) + "\n"
    assert list(er.filter_region([wide, row], "9", 1, 200))[1] == row


def test_missing_position_column_is_an_error():
    with pytest.raises(KeyError):
        list(er.filter_region(["MarkerID\tCHR\n", "a\t9\n"], "9", 0, 1))


def test_unparseable_position_is_skipped_not_fatal():
    assert kept([HEADER, "z\tchr9\tNA\t0.1\n", ROWS[2]]) == ["c"]


def test_find_merged_takes_the_largest_plain_or_gzip_table():
    base = "gs://b/sub/FelixGenome/wf/call-Concat"
    objs = [(40, f"{base}/shard-0/merged.felix.tsv.gz"), (9000, f"{base}/shard-0/attempt-2/merged.felix.tsv.gz"),
            (50, f"{base}/shard-0/stdout"), (8000, f"{base}/shard-0/other.tsv.gz")]
    assert er.find_merged(objs) == f"{base}/shard-0/attempt-2/merged.felix.tsv.gz"  # runs since 712c4a2 write .gz
    assert er.find_merged([(7, f"{base}/merged.felix.tsv")]) == f"{base}/merged.felix.tsv"  # older runs wrote the plain table
    assert er.find_merged([(7, f"{base}/stdout"), (8, f"{base}/old_merged.felix.tsv.gz")]) is None
    assert er.find_merged([]) is None


class _FakePopen:
    def __init__(self, payload: bytes):
        self.stdout = io.BytesIO(payload)

    def wait(self):
        return 0


def test_split_windows_routes_each_line_to_its_window():
    got = list(er.split_windows([HEADER] + ROWS, [("9", 8_480_000, 8_500_000), ("10", 1, 20_000_000)]))
    assert got[0][0] == -1  # the header comes first
    assert [(i, line.split("\t")[0]) for i, line in got[1:]] == [(0, "b"), (1, "g")]


def test_split_windows_yields_a_line_once_per_overlapping_window():
    got = list(er.split_windows([HEADER] + ROWS, [("9", 8_900_000, 9_000_000), ("9", 8_990_000, 9_100_000)]))
    assert [(i, line.split("\t")[0]) for i, line in got[1:]] == [(0, "c"), (1, "c")]


def test_split_windows_needs_the_position_columns():
    with pytest.raises(KeyError, match="POS"):
        list(er.split_windows(["MarkerID\tCHR\n", "a\tchr9\n"], [("9", 1, 2)]))


def test_extract_many_cuts_several_windows_in_one_pass(tmp_path, monkeypatch):
    monkeypatch.setattr(er.subprocess, "Popen", lambda *a, **k: _FakePopen((HEADER + "".join(ROWS)).encode()))
    a, b = tmp_path / "a.tsv", tmp_path / "b.tsv"
    counts = er.extract_many("gs://b/merged.felix.tsv",
                             [("9", 8_480_000, 8_500_000, a), ("9", 9_000_000, 9_480_000, b)])
    assert counts == [1, 1]
    for path, marker in ((a, "b"), (b, "d")):
        lines = path.read_text().splitlines()
        assert lines[0].startswith("MarkerID")  # every output carries the header
        assert [x.split("\t")[0] for x in lines[1:]] == [marker]


def test_extract_many_writes_a_header_only_file_for_an_empty_window(tmp_path, monkeypatch):
    monkeypatch.setattr(er.subprocess, "Popen", lambda *a, **k: _FakePopen((HEADER + "".join(ROWS)).encode()))
    out = tmp_path / "empty.tsv"
    assert er.extract_many("gs://b/merged.felix.tsv", [("22", 1, 1000, out)]) == [0]
    assert out.read_text() == HEADER


def test_read_windows_skips_header_and_comments(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text("# a comment\nchrom\tstart\tend\tout_prefix\n9\t1\t100\tloc_a\n10\t5\t50\tloc_b\n\n")
    assert er.read_windows(p) == [("9", 1, 100, "loc_a"), ("10", 5, 50, "loc_b")]


def test_read_windows_rejects_a_short_row(tmp_path):
    p = tmp_path / "w.tsv"
    p.write_text("9\t1\t100\n")
    with pytest.raises(ValueError, match="4 tab-separated fields"):
        er.read_windows(p)


@pytest.mark.parametrize("name", ["merged.felix.tsv", "merged.felix.tsv.gz"])
def test_extract_streams_plain_and_gzip_tables(tmp_path, monkeypatch, name):
    raw = (HEADER + "".join(ROWS)).encode()
    payload = gzip.compress(raw) if name.endswith(".gz") else raw
    monkeypatch.setattr(er.subprocess, "Popen", lambda *a, **k: _FakePopen(payload))
    out = tmp_path / "x.tsv"
    n = er.extract(f"gs://b/{name}", out, "9", 8_480_000, 9_480_000)
    assert n == 3
    assert [x.split("\t")[0] for x in out.read_text().splitlines()] == ["MarkerID", "b", "c", "d"]
