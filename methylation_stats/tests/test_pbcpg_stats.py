#!/usr/bin/env python3
"""Unit tests for pb-CpG-tools bedMethyl summaries."""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path

import pytest

REPO = Path(__file__).resolve().parents[2]
SCRIPTS = REPO / "scripts"
sys.path.insert(0, str(SCRIPTS))

from pbcpg_stats import (  # noqa: E402
    attr_to_str,
    combine_sample_stats,
    group_mean_concordance,
    inventory_from_frame,
    manuscript_numbers,
    paired_concordance,
    pick_concordance_samples,
    sites_uri_column,
    stats_uri_column,
    summarize_bed,
    write_concordance_entity_tsv,
)


HEADER = """\
##pb-cpg-tools-version=3.0.0
##pileup-mode=model
##modsites-mode=denovo
##min-coverage=4
##basemod-source=jasmine 2.2.0
#chrom	begin	end	mod_score	type	cov	est_mod_count	est_unmod_count	discretized_mod_score
"""


def write_bed(path: Path, rows: list[str], *, gzipped: bool = True) -> Path:
    text = HEADER + "".join(row if row.endswith("\n") else row + "\n" for row in rows)
    if gzipped:
        path.write_bytes(gzip.compress(text.encode()))
    else:
        path.write_text(text)
    return path


@pytest.fixture
def beds(tmp_path: Path) -> dict[str, Path]:
    combined = write_bed(
        tmp_path / "s1.combined.bed.gz",
        [
            "chr1\t100\t102\t90.0\tcombined\t12\t10.8\t1.2\t90",
            "chr1\t200\t202\t10.0\tcombined\t6\t0.6\t5.4\t10",
            "chr22\t50\t52\t50.0\tcombined\t4\t2.0\t2.0\t50",
            "chrM\t10\t12\t80.0\tcombined\t20\t16.0\t4.0\t80",
            "chrX\t5\t7\t70.0\tcombined\t8\t5.6\t2.4\t70",
        ],
    )
    hap1 = write_bed(
        tmp_path / "s1.hap1.bed.gz",
        [
            "chr1\t100\t102\t95.0\thap1\t6\t5.7\t0.3\t95",
            "chr22\t50\t52\t20.0\thap1\t3\t0.6\t2.4\t20",
        ],
    )
    hap2 = write_bed(
        tmp_path / "s1.hap2.bed.gz",
        [
            "chr1\t100\t102\t85.0\thap2\t6\t5.1\t0.9\t85",
            "chr1\t200\t202\t5.0\thap2\t5\t0.25\t4.75\t5",
        ],
    )
    return {"combined": combined, "hap1": hap1, "hap2": hap2}


def test_summarize_drops_mito_and_counts_thresholds(beds: dict[str, Path]) -> None:
    stats = summarize_bed(beds["combined"])
    assert stats["n_sites"] == 4  # mito dropped
    assert stats["n_autosome"] == 3
    assert stats["n_chrX"] == 1
    assert stats["n_cov_ge5"] == 3
    assert stats["n_cov_ge10"] == 1
    assert stats["mean_cov"] == pytest.approx((12 + 6 + 4 + 8) / 4)
    assert stats["pileup_mode"] == "model"
    assert "jasmine" in stats["basemod_source"]


def test_combine_haplotype_flags(beds: dict[str, Path]) -> None:
    row = combine_sample_stats(
        "1000234",
        summarize_bed(beds["combined"]),
        summarize_bed(beds["hap1"]),
        summarize_bed(beds["hap2"]),
    )
    assert row["sample_id"] == "1000234"
    assert row["haplotype_tracks"] is True
    assert row["haplotype_resolved"] is True
    assert row["combined_n_sites"] == 4
    assert row["hap1_n_sites"] == 2
    assert row["hap2_n_sites"] == 2
    assert 0 < row["frac_hap_slots"] < 1


def test_cli_summarize_and_merge(tmp_path: Path, beds: dict[str, Path]) -> None:
    from pbcpg_stats import main

    prefix = tmp_path / "shards" / "1000234"
    prefix.parent.mkdir()
    main(
        [
            "summarize",
            "--sample-id",
            "1000234",
            "--combined",
            str(beds["combined"]),
            "--hap1",
            str(beds["hap1"]),
            "--hap2",
            str(beds["hap2"]),
            "--out-prefix",
            str(prefix),
        ]
    )
    tsv = Path(str(prefix) + ".tsv")
    assert tsv.is_file()
    out_dir = tmp_path / "out"
    main(["merge", "--stats-dir", str(prefix.parent), "--out-dir", str(out_dir)])
    payload = json.loads((out_dir / "methylation_manuscript_numbers.json").read_text())
    assert payload["n_with_combined"] == 1
    assert payload["n_haplotype_resolved"] == 1
    assert payload["median_n_cpg_ge5"] == 3
    assert "Suggested" not in json.dumps(payload)


def test_paired_concordance(tmp_path: Path) -> None:
    a = write_bed(
        tmp_path / "a.bed.gz",
        [
            "chr22\t10\t12\t80.0\tcombined\t12\t9.6\t2.4\t80",
            "chr22\t20\t22\t20.0\tcombined\t12\t2.4\t9.6\t20",
        ],
    )
    b = write_bed(
        tmp_path / "b.bed.gz",
        [
            "chr22\t10\t12\t82.0\tcombined\t11\t9.0\t2.0\t82",
            "chr22\t20\t22\t50.0\tcombined\t11\t5.5\t5.5\t50",
        ],
        gzipped=True,
    )
    result = paired_concordance(a, b, min_cov=10)
    assert result["n_sites"] == 2
    assert result["concordance_frac"] == pytest.approx(0.5)


def test_group_mean_concordance(tmp_path: Path) -> None:
    p1 = tmp_path / "p1.tsv.gz"
    j1 = tmp_path / "j1.tsv.gz"
    with gzip.open(p1, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tbegin\tcov\tmod_score\nchr22\t1\t12\t80\nchr22\t2\t12\t10\n")
    with gzip.open(j1, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tbegin\tcov\tmod_score\nchr22\t1\t12\t81\nchr22\t2\t12\t40\n")
    result = group_mean_concordance(
        {"p1": p1, "j1": j1},
        {"p1": "primrose-1.4.0", "j1": "jasmine-2.2.0"},
        min_cov=10,
        min_n_per_group=1,
    )
    assert result["n_sites"] == 2
    assert result["concordance_frac"] == pytest.approx(0.5)


def test_inventory_counts_gs_uris() -> None:
    import pandas as pd

    df = pd.DataFrame(
        {
            "entity:aou2_v1_phased_bams_id": ["1", "2"],
            "combined_bed": ["gs://b/a.bed.gz", "gs://b/b.bed.gz"],
            "hap1_bed": ["gs://b/a.hap1.bed.gz", ""],
            "hap2_bed": ["gs://b/a.hap2.bed.gz", "gs://b/b.hap2.bed.gz"],
        }
    )
    result = inventory_from_frame(df, id_column="entity:aou2_v1_phased_bams_id")
    assert result["n_rows"] == 2
    assert result["n_combined_bed"] == 2
    assert result["n_haplotype_pair"] == 1
    assert result["n_stats_tsv"] == 0


def test_attr_to_str_unwraps_file_objects() -> None:
    assert attr_to_str(None) == ""
    assert attr_to_str("gs://b/a.tsv") == "gs://b/a.tsv"
    assert attr_to_str({"value": "gs://b/a.tsv"}) == "gs://b/a.tsv"
    assert attr_to_str({"path": "gs://b/a.tsv"}) == "gs://b/a.tsv"
    assert attr_to_str({"entityName": "1000234", "entityType": "aou2_v1_phased_bams"}) == "1000234"


def test_stats_uri_column_prefers_stats_tsv() -> None:
    import pandas as pd

    df = pd.DataFrame(
        {
            "entity:aou2_v1_phased_bams_id": ["1"],
            "combined_bed": ["gs://b/a.bed.gz"],
            "stats_tsv": ["gs://b/1.stats.tsv"],
        }
    )
    assert stats_uri_column(df) == "stats_tsv"
    assert inventory_from_frame(df, id_column="entity:aou2_v1_phased_bams_id")["n_stats_tsv"] == 1


def test_pick_concordance_samples_balances_callers(tmp_path: Path) -> None:
    import pandas as pd

    df = pd.DataFrame(
        {
            "entity:aou2_v1_phased_bams_id": [str(i) for i in range(10)],
            "combined_bed": [f"gs://b/{i}.bed.gz" for i in range(10)],
        }
    )
    cov = tmp_path / "cov.csv"
    pd.DataFrame(
        {
            "research_id": [str(i) for i in range(10)],
            "pb_meth_caller": ["primrose-1.4.0"] * 6 + ["jasmine-2.2.0"] * 3 + ["mixed"],
        }
    ).to_csv(cov, index=False)
    stats = tmp_path / "stats.tsv"
    pd.DataFrame(
        {
            "sample_id": [str(i) for i in range(10)],
            "combined_mean_cov": list(range(10)),
        }
    ).to_csv(stats, sep="\t", index=False)
    picked = pick_concordance_samples(
        df,
        covariates=cov,
        stats=stats,
        n_per_group=2,
    )
    assert list(picked["caller_group"]).count("primrose") == 2
    assert list(picked["caller_group"]).count("jasmine") == 2
    assert "mixed" not in set(picked["caller_group"])
    primrose_ids = set(
        picked.loc[picked["caller_group"] == "primrose", "entity:aou2_v1_phased_bams_id"]
    )
    jasmine_ids = set(
        picked.loc[picked["caller_group"] == "jasmine", "entity:aou2_v1_phased_bams_id"]
    )
    assert primrose_ids == {"4", "5"}
    assert jasmine_ids == {"7", "8"}
    entity_path = tmp_path / "entities.tsv"
    write_concordance_entity_tsv(
        picked, entity_path, id_column="entity:aou2_v1_phased_bams_id"
    )
    header = entity_path.read_text().splitlines()[0]
    assert header.startswith("entity:pbcpg_concordance_id")


def test_sites_uri_column() -> None:
    import pandas as pd

    df = pd.DataFrame(
        {
            "entity:aou2_v1_phased_bams_id": ["1"],
            "combined_bed": ["gs://b/a.bed.gz"],
            "sites": ["gs://b/1.sites.tsv.gz"],
        }
    )
    assert sites_uri_column(df) == "sites"
    inv = inventory_from_frame(df, id_column="entity:aou2_v1_phased_bams_id")
    assert inv["n_sites"] == 1
    assert inv["sites_uri_column"] == "sites"


def test_entities_to_rows_unwraps_file_attrs() -> None:
    from pbcpg_stats import _entities_to_rows

    rows = _entities_to_rows(
        [
            {
                "name": "1000234",
                "attributes": {
                    "combined_bed": {"path": "gs://b/a.combined.bed.gz"},
                    "hap1_bed": {"value": "gs://b/a.hap1.bed.gz"},
                    "hap2_bed": "gs://b/a.hap2.bed.gz",
                },
            }
        ],
        id_column="entity:aou2_v1_phased_bams_id",
    )
    assert rows[0]["entity:aou2_v1_phased_bams_id"] == "1000234"
    assert rows[0]["combined_bed"] == "gs://b/a.combined.bed.gz"
    assert rows[0]["hap1_bed"] == "gs://b/a.hap1.bed.gz"


def test_fetch_phased_bams_table_prefers_entities_tsv(monkeypatch) -> None:
    import sys
    import types

    import pbcpg_stats

    fake_api = types.ModuleType("firecloud.api")

    def get_entities_tsv(namespace, workspace, etype, model="flexible"):
        assert namespace == "allofus-drc-wgs-LR-prodData"
        assert workspace == "AoU_DRC_LongReads_PhaseTwo_Storage"
        assert etype == "aou2_v1_phased_bams"
        assert model == "flexible"
        return types.SimpleNamespace(
            status_code=200,
            text=(
                "entity:aou2_v1_phased_bams_id\tcombined_bed\thap1_bed\thap2_bed\n"
                "1000234\tgs://b/a.combined.bed.gz\tgs://b/a.hap1.bed.gz\tgs://b/a.hap2.bed.gz\n"
            ),
        )

    def get_entities(*args, **kwargs):
        raise AssertionError("JSON get_entities should not run when TSV succeeds")

    fake_api.get_entities_tsv = get_entities_tsv
    fake_api.get_entities = get_entities
    fake_pkg = types.ModuleType("firecloud")
    fake_pkg.api = fake_api
    monkeypatch.setitem(sys.modules, "firecloud", fake_pkg)
    monkeypatch.setitem(sys.modules, "firecloud.api", fake_api)

    df = pbcpg_stats.fetch_phased_bams_table(from_firecloud=True)
    assert list(df["entity:aou2_v1_phased_bams_id"]) == ["1000234"]
    assert df.loc[0, "combined_bed"] == "gs://b/a.combined.bed.gz"


def test_extract_chrom_window(tmp_path: Path, beds: dict[str, Path]) -> None:
    from pbcpg_stats import main

    out = tmp_path / "chr1.window.tsv"
    main(
        [
            "extract-chrom",
            "--bed",
            str(beds["combined"]),
            "--chrom",
            "chr1",
            "--start",
            "150",
            "--end",
            "250",
            "--min-cov",
            "1",
            "--out",
            str(out),
        ]
    )
    lines = [ln for ln in out.read_text().splitlines() if not ln.startswith("#") and not ln.startswith("chrom")]
    assert len(lines) == 1
    assert lines[0].startswith("chr1\t200\t")


def test_concordance_writes_pairs(tmp_path: Path) -> None:
    from pbcpg_stats import main

    p1 = tmp_path / "p1.tsv.gz"
    j1 = tmp_path / "j1.tsv.gz"
    with gzip.open(p1, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tbegin\tcov\tmod_score\nchr22\t1\t12\t80\nchr22\t2\t12\t10\n")
    with gzip.open(j1, "wt", encoding="utf-8") as handle:
        handle.write("chrom\tbegin\tcov\tmod_score\nchr22\t1\t12\t81\nchr22\t2\t12\t40\n")
    labels = tmp_path / "labels.tsv"
    labels.write_text("sample_id\tpb_meth_caller\np1\tprimrose-1.4.0\nj1\tjasmine-2.2.0\n")
    out = tmp_path / "methylation_concordance.json"
    main(
        [
            "concordance",
            "--dumps-dir",
            str(tmp_path),
            "--labels",
            str(labels),
            "--min-n-per-group",
            "1",
            "--out",
            str(out),
        ]
    )
    pairs = tmp_path / "methylation_concordance.pairs.tsv.gz"
    assert pairs.is_file()
    payload = json.loads(out.read_text())
    assert payload["n_sites"] == 2
    assert payload["concordance_frac"] == pytest.approx(0.5)


def test_plot_s3(tmp_path: Path) -> None:
    pytest.importorskip("matplotlib")
    from pbcpg_stats import main

    pairs = tmp_path / "pairs.tsv"
    pairs.write_text("mean_primrose\tmean_jasmine\n80\t81\n10\t12\n50\t48\n")
    stats = tmp_path / "pbcpg_sample_stats.tsv"
    stats.write_text(
        "sample_id\tcoverage\tfrac_hap_slots_ge5\n"
        "1\t16.0\t0.70\n"
        "2\t18.0\t0.74\n"
        "3\t32.0\t0.88\n"
    )
    prefix = tmp_path / "figS3"
    main(
        [
            "plot-s3",
            "--pairs",
            str(pairs),
            "--stats",
            str(stats),
            "--out-prefix",
            str(prefix),
        ]
    )
    assert prefix.with_suffix(".png").is_file() or Path(str(prefix) + ".png").is_file()
    assert Path(str(prefix) + ".png").is_file()
    assert Path(str(prefix) + ".pdf").is_file()
