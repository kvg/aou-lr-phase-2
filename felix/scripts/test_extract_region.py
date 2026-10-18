"""Tests for the window filter in extract_region.py (the gsutil calls are not exercised here)."""

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
