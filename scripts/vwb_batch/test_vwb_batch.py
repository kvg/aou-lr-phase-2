#!/usr/bin/env python3
"""Local checks for vwb_batch (no GCP)."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from vwb_batch import alloc, cohort, cost, monitor, plots, submit, wave  # noqa: E402
from vwb_batch.pipelines import expansion_hunter as eh  # noqa: E402
from vwb_batch.pipelines import locityper as lt  # noqa: E402


def _records(n: int) -> list[dict[str, str]]:
    rows = []
    for i in range(n):
        sid = f"s{i:06d}"
        rows.append(
            {
                "sample_id": sid,
                "cram": f"gs://bucket/wgs_{sid}.cram",
                "crai": f"gs://bucket/wgs_{sid}.cram.crai",
                "ref_fa": "gs://bucket/ref.fa",
                "ref_fai": "gs://bucket/ref.fa.fai",
                "catalog": "gs://bucket/catalog.json",
                "sex": "female",
            }
        )
    return rows


def test_cost() -> None:
    xs = [1.0, 2.0, 3.0, 4.0, 100.0]
    assert abs(cost.quantile(xs, 0.5) - 3.0) < 1e-9
    p90 = cost.quantile(xs, 0.9)
    assert p90 > 4.0
    parsed = cost.parse_transfer_stats(
        "file_path\ttotal_bytes_loaded\ttotal_containers_loaded_from_cram\t"
        "total_read_requests\ttotal_index_bytes_loaded\ttotal_network_bytes\n"
        "gs://x\t10\t1\t1\t2\t1073741824\n"
    )
    assert parsed == 1073741824
    near = cost.nearline_cost_usd(parsed)
    assert abs(near - 0.01) < 1e-9
    sample = cost.sample_cost_usd(
        minicram_seconds=3600, genotype_seconds=1800, network_bytes=parsed
    )
    assert sample["total_usd"] > 0.15
    summary = cost.summarize_pilot([sample, sample], q=0.9)
    n = cost.n_keep(budget_usd=100.0, planning_usd=summary["planning_usd"], retry_margin=1.15, n_available=10_000)
    assert 0 < n < 10_000
    print("cost ok", round(sample["total_usd"], 4), "n_keep", n)


def test_job_json_tsv() -> None:
    records = _records(500)
    rows = eh.env_rows(records, stage="minicram", out_prefix="gs://out/eh", project="proj")
    job = eh.build_job(
        stage="minicram",
        rows=rows,
        tasks_tsv_uri="gs://out/runs/r/shards/minicram/shard-0000/tasks.tsv",
        out_prefix="gs://out/eh",
        project="proj",
        region="us-central1",
        sa="pet@example.com",
        image=eh.PRINT_READS_DOCKER,
        parallelism=50,
        labels=alloc.job_labels(pipeline="expansion_hunter", run_id="eh-test", stage="minicram", shard=0),
        max_run_duration="28800s",
    )
    group = job["taskGroups"][0]
    assert group["taskCount"] == 500
    assert group["parallelism"] == 50
    assert "taskEnvironments" not in group
    assert group["taskSpec"]["environment"]["variables"]["TASKS_TSV"].endswith("tasks.tsv")
    size = submit.assert_job_json_ok(job)
    assert size < submit.MAX_JOB_JSON_BYTES
    assert job["allocationPolicy"]["network"]["networkInterfaces"][0]["noExternalIpAddress"] is True
    print("minicram TASKS_TSV job json", size, "bytes")

    gt_rows = eh.env_rows(records[:80], stage="genotype", out_prefix="gs://out/eh", project="proj")
    gt = eh.build_job(
        stage="genotype",
        rows=gt_rows,
        tasks_tsv_uri="gs://out/runs/r/shards/genotype/shard-0000/tasks.tsv",
        out_prefix="gs://out/eh",
        project="proj",
        region="us-central1",
        sa="pet@example.com",
        image=eh.EH_DOCKER,
        parallelism=20,
        labels=alloc.job_labels(pipeline="expansion_hunter", run_id="eh-test", stage="genotype", shard=0),
    )
    assert "taskEnvironments" not in gt["taskGroups"][0]
    assert gt["taskGroups"][0]["taskSpec"]["runnables"][2]["alwaysRun"] is True
    print("genotype TASKS_TSV job json", submit.job_json_size(gt), "bytes")


def test_smoke_csv_job_still_embeds_env() -> None:
    csv_path = ROOT / "expansion_hunter" / "configs" / "batch.header.csv"
    mod = eh._load_submit_batch()
    job = mod.build_minicram_job(
        csv_path=csv_path,
        worker_path=ROOT / "expansion_hunter" / "batch" / "make_minicram.py",
        image=eh.PRINT_READS_DOCKER,
        out_prefix="gs://out/eh",
        project="proj",
        region="us-central1",
        sa="pet@example.com",
    )
    env = job["taskGroups"][0]["taskSpec"]["environment"]["variables"]
    assert env["SAMPLE_ID"] == "1000000"
    assert env["CRAM"].startswith("gs://")
    print("smoke csv embed ok")


def test_shards_and_ids() -> None:
    parts = cohort.shards(_records(250), 100)
    assert [len(p[1]) for p in parts] == [100, 100, 50]
    jid = alloc.job_id(pipeline="expansion_hunter", stage="minicram", run_id="eh-260929-pilot", shard=3)
    assert jid.startswith("eh-mc-")
    assert len(jid) <= 63
    with tempfile.TemporaryDirectory() as tmp:
        path = Path(tmp) / "t.tsv"
        cohort.write_tasks_tsv(path, [{"SAMPLE_ID": "a", "CRAM": "gs://x"}])
        assert "SAMPLE_ID" in path.read_text()
    print("shards/job_id ok", jid)


def test_recommend_compute() -> None:
    rows = [
        {
            "peak_rss_bytes": 1.5 * 1024**3,
            "peak_disk_bytes": 4 * 1024**3,
            "cpu_cores_avg": 0.4,
            "wall_sec": 400,
            "cpu_sec": 160,
        }
        for _ in range(10)
    ]
    rec = cost.recommend_compute(rows, stage="minicram")
    assert rec["cpuMilli"] == 1000
    assert rec["memoryMib"] == 3072  # 1.5 GiB * 1.4 → 2.1 → ceil to 3 GiB
    assert rec["bootDiskMib"] == 20480
    assert rec["machine"] == "e2-standard-2"
    fat = cost.recommend_compute(
        [
            {
                "peak_rss_bytes": 20 * 1024**3,
                "peak_disk_bytes": 40 * 1024**3,
                "cpu_cores_avg": 3.5,
                "wall_sec": 1,
                "cpu_sec": 3.5,
            }
        ],
        stage="minicram",
    )
    assert fat["memoryMib"] == 16384  # cap at default
    assert fat["cpuMilli"] == 4000
    print("recommend_compute ok", rec)


def test_worker_source_compiles() -> None:
    mod = eh._load_submit_batch()
    src = mod._worker_source(ROOT / "expansion_hunter" / "batch" / "make_minicram.py")
    compile(src, "minicram-worker", "exec")
    assert "class ResourceMonitor" in src
    job = mod.build_minicram_job(
        csv_path=ROOT / "expansion_hunter" / "configs" / "batch.header.csv",
        worker_path=ROOT / "expansion_hunter" / "batch" / "make_minicram.py",
        image=eh.PRINT_READS_DOCKER,
        out_prefix="gs://out/eh",
        project="proj",
        region="us-central1",
        sa="pet@example.com",
        compute={"cpuMilli": 1000, "memoryMib": 4096, "bootDiskMib": 20480},
    )
    cr = job["taskGroups"][0]["taskSpec"]["computeResource"]
    assert cr == {"cpuMilli": 1000, "memoryMib": 4096, "bootDiskMib": 20480}
    env = job["taskGroups"][0]["taskSpec"]["environment"]["variables"]
    assert env["RESOURCE_STATS"].endswith("minicram.resources.tsv")
    print("worker source + tuned compute ok")


def test_locityper_job_json() -> None:
    records = []
    for i in range(80):
        sid = f"s{i:06d}"
        records.append(
            {
                "sample_id": sid,
                "cram": f"gs://bucket/wgs_{sid}.cram",
                "crai": f"gs://bucket/wgs_{sid}.cram.crai",
                "ref_fa": "gs://bucket/ref.fa",
                "ref_fai": "gs://bucket/ref.fa.fai",
                "counts_jf": "gs://bucket/counts.jf",
                "bed": "gs://bucket/loci.bed",
                "db_tar": "gs://bucket/vcf_db.tar.gz",
                "sex": "female",
            }
        )
    rows = lt.env_rows(records, stage="minicram", out_prefix="gs://out/lt", project="proj")
    job = lt.build_job(
        stage="minicram",
        rows=rows,
        tasks_tsv_uri="gs://out/runs/r/shards/minicram/shard-0000/tasks.tsv",
        out_prefix="gs://out/lt",
        project="proj",
        region="us-central1",
        sa="pet@example.com",
        image=lt.PRINT_READS_DOCKER,
        parallelism=20,
        labels=alloc.job_labels(pipeline="locityper", run_id="lt-test", stage="minicram", shard=0),
    )
    group = job["taskGroups"][0]
    assert group["taskCount"] == 80
    assert "taskEnvironments" not in group
    size = submit.assert_job_json_ok(job)
    assert size < submit.MAX_JOB_JSON_BYTES
    print("locityper minicram TASKS_TSV job json", size, "bytes")

    gt_rows = lt.env_rows(records[:10], stage="genotype", out_prefix="gs://out/lt", project="proj")
    gt = lt.build_job(
        stage="genotype",
        rows=gt_rows,
        tasks_tsv_uri="gs://out/runs/r/shards/genotype/shard-0000/tasks.tsv",
        out_prefix="gs://out/lt",
        project="proj",
        region="us-central1",
        sa="pet@example.com",
        image=lt.LOCITYPER_DOCKER,
        parallelism=10,
        labels=alloc.job_labels(pipeline="locityper", run_id="lt-test", stage="genotype", shard=0),
    )
    spec = gt["taskGroups"][0]["taskSpec"]
    ctr = spec["runnables"][1]["container"]
    assert ctr["entrypoint"] == "/bin/bash"
    assert "locityper preproc" in ctr["commands"][1]
    assert "python3" in spec["runnables"][2]["script"]["text"]
    assert "LT_WORKDIR" in spec["runnables"][2]["script"]["text"]
    assert spec["runnables"][3]["alwaysRun"] is True
    compile((ROOT / "locityper" / "batch" / "genotype.py").read_text(), "lt-genotype-summary", "exec")
    assert spec["computeResource"]["cpuMilli"] == 2000
    csv_path = ROOT / "locityper" / "configs" / "batch.header.csv"
    mod = lt._load_submit_batch()
    smoke = mod.build_minicram_job(
        csv_path=csv_path,
        worker_path=ROOT / "locityper" / "batch" / "make_minicram.py",
        image=lt.PRINT_READS_DOCKER,
        out_prefix="gs://out/lt",
        project="proj",
        region="us-central1",
        sa="pet@example.com",
    )
    env = smoke["taskGroups"][0]["taskSpec"]["environment"]["variables"]
    assert env["SAMPLE_ID"] == "1000000"
    assert env["BED"].endswith("smoke.bed")
    src = mod._worker_source(ROOT / "locityper" / "batch" / "make_minicram.py")
    compile(src, "lt-minicram-worker", "exec")
    assert "HOST_LOG_GCS" in gt_rows[0]
    assert "host.log" in mod.UPLOAD_SCRIPT
    assert "pkg.dev/" in lt.LOCITYPER_DOCKER
    try:
        mod.require_artifact_registry("eichlerlab/locityper:1.4.5.0")
    except SystemExit:
        pass
    else:
        raise AssertionError("docker hub image should be rejected")
    mod.require_artifact_registry(lt.PRINT_READS_DOCKER)
    print("locityper genotype + smoke csv ok", submit.job_json_size(gt), "bytes")


def test_monitor_summary() -> None:
    job = {
        "name": "projects/p/locations/us-central1/jobs/eh-mc-pilot-s000-120000",
        "createTime": "2026-09-28T20:00:00Z",
        "labels": {
            "pipeline": "expansion-hunter",
            "run-id": "eh-260928-pilot",
            "stage": "minicram",
            "shard": "0",
        },
        "taskGroups": [
            {
                "taskCount": 20,
                "parallelism": 20,
                "taskSpec": {"computeResource": {"cpuMilli": 4000, "memoryMib": 16384}},
            }
        ],
        "allocationPolicy": {"instances": [{"policy": {"machineType": "e2-standard-4"}}]},
        "status": {
            "state": "RUNNING",
            "runDuration": "120s",
            "taskGroups": {"group0": {"counts": {"SUCCEEDED": "5", "RUNNING": "15"}}},
        },
    }
    row = monitor.summarize_job(job)
    assert row["job_id"] == "eh-mc-pilot-s000-120000"
    assert row["succeeded"] == 5
    assert row["running"] == 15
    assert row["inflight_vcpu"] == 60
    assert row["age"] == "2m00s"
    done = dict(job)
    done["status"] = {
        "state": "SUCCEEDED",
        "runDuration": "480s",
        "taskGroups": {"group0": {"counts": {"SUCCEEDED": "20"}}},
    }
    assert monitor.summarize_job(done)["inflight_vcpu"] == 0
    kept = monitor.filter_jobs([job], name_prefix="eh-", run_id="eh-260928-pilot")
    assert len(kept) == 1
    assert monitor.filter_jobs([job], name_prefix="lt-") == []
    table = monitor.format_table([row])
    assert "eh-mc-pilot-s000-120000" in table
    print("monitor ok", row["inflight_vcpu"], "vCPU")


def test_wave_prefers_first_attempts() -> None:
    ids = ["new", "stranded", "failed", "running", "done", "missing"]
    events = [
        {
            "event": "job_submitted",
            "stage": "minicram",
            "job_id": "job-a",
            "sample_ids": ["stranded", "failed"],
        },
        {
            "event": "job_submitted",
            "stage": "minicram",
            "job_id": "job-b",
            "sample_ids": ["running", "done"],
        },
        {
            "event": "job_submitted",
            "stage": "minicram",
            "job_id": "job-c",
            "sample_ids": ["missing"],
        },
    ]
    polls = {
        "job-a": {
            "state": "FAILED",
            "tasks": [
                {"index": 0, "state": "PENDING"},
                {"index": 1, "state": "FAILED"},
            ],
        },
        "job-b": {
            "state": "RUNNING",
            "tasks": [
                {"index": 0, "state": "RUNNING"},
                {"index": 1, "state": "SUCCEEDED"},
            ],
        },
        "job-c": {
            "state": "SUCCEEDED",
            "tasks": [{"index": 0, "state": "SUCCEEDED"}],
        },
    }
    cats = wave.classify_stage(
        ids,
        submissions=wave.latest_submissions(events, "minicram"),
        polls=polls,
        done_ids={"done"},
    )
    assert cats["new"] == "not_started"
    assert cats["stranded"] == "not_started"
    assert cats["failed"] == "failed"
    assert cats["running"] == "running"
    assert cats["done"] == "done"
    assert cats["missing"] == "failed"
    assert wave.choose_to_submit(ids, cats) == ["new", "stranded"]
    rest = wave.choose_to_submit(["failed", "missing", "done"], cats)
    assert rest == ["failed", "missing"]
    assert wave.count_categories(cats)["not_started"] == 2
    print("wave ok", wave.count_categories(cats))


def test_progress_figure() -> None:
    try:
        import matplotlib.pyplot as plt
    except ImportError:
        print("progress figure skipped (no matplotlib)")
        return
    counts = {
        "minicram": {"done": 20, "running": 800, "failed": 3, "not_started": 1012},
        "genotype": {"done": 20, "running": 0, "failed": 0, "not_started": 1815},
    }
    with tempfile.TemporaryDirectory() as tmp:
        paths = plots.plot_cohort_progress(
            Path(tmp),
            run_id="eh-test",
            catalog="candidate_EH_Loci.GRCh38.degenerate.json",
            counts=counts,
            fully_done=20,
            keep_n=1835,
            inflight_vcpu=800,
            e2_quota=2400,
            n_active_jobs=4,
            updated_at="2026-09-29T00:00:00Z",
        )
        assert len(paths) == 2
        assert all(p.is_file() and p.stat().st_size > 0 for p in paths)
    plt.close("all")
    print("progress figure ok")


def main() -> None:
    test_cost()
    test_job_json_tsv()
    test_smoke_csv_job_still_embeds_env()
    test_shards_and_ids()
    test_recommend_compute()
    test_worker_source_compiles()
    test_locityper_job_json()
    test_monitor_summary()
    test_wave_prefers_first_attempts()
    test_progress_figure()
    print("all ok")


if __name__ == "__main__":
    main()
