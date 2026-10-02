#!/usr/bin/env python3
"""Local checks for vwb_batch (no GCP)."""

from __future__ import annotations

import json
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "scripts"))

from vwb_batch import alloc, campaign, cohort, cost, monitor, plots, submit, wave  # noqa: E402
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


def test_fused_job_json() -> None:
    recs = _records(1000)
    rows = eh.env_rows(recs, stage="fused", out_prefix="gs://out/v2", project="proj")
    assert rows[0]["EH_VCF"] == "gs://out/v2/s000000/s000000.EH.vcf"
    assert rows[0]["MINICRAM"].endswith("s000000.minicram.cram")
    job = eh.build_fused_job(
        rows=rows,
        tasks_tsv_uri="gs://ledger/tasks/w001/j.tsv",
        project="proj",
        region="us-central1",
        sa="pet@proj.iam.gserviceaccount.com",
        labels={"run-id": "r"},
        parallelism=1000,
        spot=True,
    )
    size = submit.assert_job_json_ok(job)
    group = job["taskGroups"][0]
    spec = group["taskSpec"]
    assert group["taskCount"] == 1000 and "taskEnvironments" not in group
    assert "taskCountPerNode" not in group  # packing is safe: per-task directories
    assert spec["environment"]["variables"]["EH_TASK_ROOT"] == "/mnt/disks/eh"
    assert spec["environment"]["variables"]["TASKS_TSV"] == "gs://ledger/tasks/w001/j.tsv"
    run = spec["runnables"]
    assert [("script" in r, "container" in r) for r in run] == [(True, False), (False, True), (False, True), (True, False)]
    assert run[-1].get("alwaysRun") is True
    assert all(r["container"]["volumes"] == ["/mnt/disks/eh:/mnt/disks/eh"] for r in run[1:3])
    for r in (run[0], run[-1]):
        assert 't${BATCH_TASK_INDEX}' in r["script"]["text"]
    assert spec["lifecyclePolicies"][0]["actionCondition"]["exitCodes"] == [50001, 50002]
    assert spec["maxRetryCount"] == 2
    assert job["allocationPolicy"]["instances"] == [{"policy": {"provisioningModel": "SPOT"}}]
    assert spec["computeResource"] == {"cpuMilli": 1000, "memoryMib": 2048, "bootDiskMib": 20480}
    for r in run[1:3]:
        compile(r["container"]["commands"][1], "<worker>", "exec")
    legacy = eh.build_job(
        stage="genotype",
        rows=eh.env_rows(_records(2), stage="genotype", out_prefix="gs://o", project="p"),
        tasks_tsv_uri="gs://l/t.tsv",
        out_prefix="gs://o",
        project="p",
        region="us-central1",
        sa="sa",
        image="img",
        parallelism=2,
        labels={},
    )
    assert legacy["taskGroups"][0]["taskCountPerNode"] == 1
    assert len(eh.code_version()) == 10
    print("fused job ok", size, "bytes")


def _job(job_id, sids, *, kind="wave", wave=1, cv="v1", states=None, running=False, counts=None, cost_usd=0.0):
    entry = {
        "job_id": job_id,
        "wave": wave,
        "kind": kind,
        "code_version": cv,
        "submitted_at": f"2026-10-03T00:{wave:02d}:00Z",
        "sample_ids": list(sids),
        "state": "RUNNING" if running else "SUCCEEDED",
    }
    if not running:
        entry["task_states"] = list(states or ["SUCCEEDED"] * len(sids))
        entry["cost"] = {"vm_usd": cost_usd, "nearline_usd": 0.0, "n_ran": float(len(sids))}
    else:
        entry["counts"] = counts or {"RUNNING": len(sids)}
    return entry


def test_campaign_plan() -> None:
    ids = [f"s{i:06d}" for i in range(30_000)]
    cfg = campaign.Settings(out_prefix="gs://o", budget_usd=3000, wave_size=10_000, canary_n=200, job_size=1000)

    def state(*jobs):
        return {"jobs": {j["job_id"]: j for j in jobs}}

    p = campaign.plan(cfg, ids, state(), code_version="v1")
    assert p.action == "submit" and p.kind == "canary" and p.sample_ids == ids[:200], p.explain()

    canary_running = _job("c", ids[:200], kind="canary", running=True)
    p = campaign.plan(cfg, ids, state(canary_running), code_version="v1")
    assert p.action == "wait" and "canary" in p.reason

    bad = ["FAILED"] * 30 + ["SUCCEEDED"] * 170
    p = campaign.plan(cfg, ids, state(_job("c", ids[:200], kind="canary", states=bad)), code_version="v1")
    assert p.action == "halt" and "15%" in p.reason, p.reason

    # A code change restarts the canary even after a bad one.
    p = campaign.plan(cfg, ids, state(_job("c", ids[:200], kind="canary", states=bad)), code_version="v2")
    assert p.action == "submit" and p.kind == "canary"
    assert p.sample_ids[:30] == ids[200:230], "first attempts go before retries"

    lost_canary = _job("c", ids[:200], kind="canary", states=[""] * 200)
    p = campaign.plan(cfg, ids, state(lost_canary), code_version="v1")
    assert p.action == "submit" and p.kind == "canary" and "never ran" in p.reason, p.explain()

    canary_ok = _job("c", ids[:200], kind="canary", cost_usd=4.0)
    p = campaign.plan(cfg, ids, state(canary_ok), code_version="v1")
    assert p.action == "submit" and p.kind == "wave" and len(p.sample_ids) == 10_000
    assert p.sample_ids[0] == ids[200]

    w1 = _job("w1", ids[200:6200], wave=2, running=True)
    p = campaign.plan(cfg, ids, state(canary_ok, w1), code_version="v1")
    assert p.action == "wait" and "6,000" in p.reason

    w1 = _job("w1", ids[200:4200], wave=2, running=True)
    p = campaign.plan(cfg, ids, state(canary_ok, w1), code_version="v1")
    assert p.action == "submit" and len(p.sample_ids) == 6_000 and p.sample_ids[0] == ids[4200]

    # Running jobs with too many failures stop the next wave.
    w1 = _job("w1", ids[200:4200], wave=2, running=True, counts={"FAILED": 300, "SUCCEEDED": 700, "RUNNING": 3000})
    p = campaign.plan(cfg, ids, state(canary_ok, w1), code_version="v1")
    assert p.action == "halt", p.explain()

    # Budget caps the wave: $3000 - $4 spent - 0 in flight = 119,840 at $0.025, so try a small budget.
    small = campaign.Settings(out_prefix="gs://o", budget_usd=104.0, wave_size=10_000, canary_n=200)
    p = campaign.plan(small, ids, state(canary_ok), code_version="v1")
    assert p.action == "submit" and len(p.sample_ids) == 4000 and "budget" in p.reason, p.explain()
    tiny = campaign.Settings(out_prefix="gs://o", budget_usd=4.0, wave_size=10_000, canary_n=200)
    p = campaign.plan(tiny, ids, state(canary_ok), code_version="v1")
    assert p.action == "halt" and p.reason.startswith("budget")

    # Retries only once every sample has had a first attempt; give up after max_attempts.
    few = ids[:300]
    st = [ "FAILED" if i < 10 else "SUCCEEDED" for i in range(100)]
    done = state(
        _job("c", few[:200], kind="canary"),
        _job("w1", few[200:], wave=2, states=st),
    )
    p = campaign.plan(cfg, few, done, code_version="v1")
    assert p.action == "submit" and p.kind == "retry" and p.sample_ids == few[200:210]
    r1 = _job("r1", few[200:210], kind="retry", wave=3, states=["FAILED"] * 10)
    r2 = _job("r2", few[200:210], kind="retry", wave=4, states=["FAILED"] * 9 + ["SUCCEEDED"])
    p = campaign.plan(cfg, few, state(*done["jobs"].values(), r1, r2), code_version="v1")
    assert p.action == "finished" and "9 gave up" in p.reason, p.explain()

    # Tasks that never ran (cancelled or lost job) go back to not_started without using an attempt.
    lost = _job("w1", few[200:], wave=2, states=[""] * 100)
    lost["state"] = "LOST"
    states = campaign.sample_states(few, state(_job("c", few[:200], kind="canary"), lost), max_attempts=3)
    assert states[few[250]].status == "not_started" and states[few[250]].attempts == 0
    print("campaign plan ok")


def test_job_cost_model() -> None:
    tasks = [
        {"status": {"state": "SUCCEEDED", "statusEvents": [
            {"taskState": "RUNNING", "eventTime": "2026-09-29T07:44:03.505511122Z",
             "description": "Task state is updated from ASSIGNED to RUNNING on zones/us-central1-a/instances/111"},
            {"taskState": "SUCCEEDED", "eventTime": "2026-09-29T07:54:03.505511122Z",
             "description": "Task state is updated from RUNNING to SUCCEEDED on zones/us-central1-a/instances/111"},
        ]}, "name": "projects/p/locations/l/jobs/j/taskGroups/group0/tasks/0"},
        {"status": {"state": "FAILED", "statusEvents": [
            {"taskState": "RUNNING", "eventTime": "2026-09-29T07:44:03Z",
             "description": "Task state is updated from ASSIGNED to RUNNING on zones/us-central1-a/instances/111"},
            {"taskState": "FAILED", "eventTime": "2026-09-29T07:44:13Z",
             "description": "Task state is updated from RUNNING to FAILED on zones/us-central1-a/instances/111"},
        ]}, "name": "projects/p/locations/l/jobs/j/taskGroups/group0/tasks/1"},
    ]
    parsed = [campaign.parse_task(t) for t in tasks]
    assert parsed[0]["instance"] == "111" and parsed[1]["index"] == 1 and parsed[1]["state"] == "FAILED"
    c = cost.job_cost_usd(parsed, machine="e2-highcpu-2", spot=False)
    assert c["n_vms"] == 1 and c["n_ran"] == 2 and c["n_nearline"] == 1
    rate = cost.vm_hourly_usd("e2-highcpu-2")
    assert abs(c["vm_usd"] - (600.505511 + cost.VM_STARTUP_S) / 3600 * rate) < 1e-9
    spot = cost.job_cost_usd(parsed, machine="e2-highcpu-2", spot=True)
    assert spot["vm_usd"] < c["vm_usd"]
    assert cost.machine_shape("e2-highcpu-2") == (2, 2.0, "e2")
    assert cost.machine_shape("n2-standard-4") == (4, 16.0, "n2")
    print("job cost ok", round(c["total_usd"], 4))


def test_campaign_tick_end_to_end() -> None:
    """tick() against in-memory GCS and Batch: canary, wave, failures, recovery."""
    from vwb_batch import gcs as gcs_mod, registry

    store: dict[str, str] = {}
    batch: dict[str, dict] = {}
    saved = {
        name: getattr(gcs_mod, name)
        for name in ("exists", "cat", "upload_text", "upload_text_if_absent", "delete")
    }
    saved_submit = (submit.list_jobs, submit.list_tasks, submit.submit_job)

    def upload_if_absent(uri, text):
        if uri in store:
            return False
        store[uri] = text
        return True

    def fake_submit(*, job_id, job, project, region, quiet=False):
        assert job_id not in batch
        batch[job_id] = {"job": job, "state": "QUEUED", "n": job["taskGroups"][0]["taskCount"], "fail": set()}
        return job_id

    def fake_list_jobs(*, project, region, filter_expr=""):
        out = []
        for job_id, b in batch.items():
            n_fail = len(b["fail"]) if b["state"] in ("SUCCEEDED", "FAILED") else 0
            out.append({
                "name": f"projects/p/locations/r/jobs/{job_id}",
                "labels": b["job"]["labels"],
                "createTime": "2026-10-03T00:00:00.000Z",
                "taskGroups": b["job"]["taskGroups"],
                "status": {
                    "state": b["state"],
                    "taskGroups": {"group0": {
                        "counts": {"SUCCEEDED": str(b["n"] - n_fail), "FAILED": str(n_fail)},
                        "instances": [{"machineType": "e2-highcpu-2", "provisioningModel": "SPOT", "taskPack": "2",
                                       "bootDisk": {"sizeGb": "53"}}],
                    }},
                },
            })
        return out

    def fake_list_tasks(*, job_id, project, region, quiet=False):
        b = batch[job_id]
        rows = []
        for i in range(b["n"]):
            st = "FAILED" if i in b["fail"] else "SUCCEEDED"
            inst = f"vm{i // 2}"
            rows.append({"name": f"x/tasks/{i}", "status": {"state": st, "statusEvents": [
                {"taskState": "RUNNING", "eventTime": "2026-10-03T01:00:00Z", "description": f"on zones/z/instances/{inst}"},
                {"taskState": st, "eventTime": "2026-10-03T01:13:00Z", "description": f"on zones/z/instances/{inst}"},
            ]}})
        return rows

    gcs_mod.exists = lambda uri: uri in store
    gcs_mod.cat = lambda uri: store[uri]
    gcs_mod.upload_text = lambda uri, text: store.__setitem__(uri, text)
    gcs_mod.upload_text_if_absent = upload_if_absent
    gcs_mod.delete = lambda uri: store.pop(uri, None)
    submit.list_jobs, submit.list_tasks, submit.submit_job = fake_list_jobs, fake_list_tasks, fake_submit
    try:
        paths = registry.run_paths(output_bucket="gs://b", pipeline="expansion_hunter", run_id="eh-test")
        cfg = campaign.Settings(out_prefix="gs://b/out", budget_usd=3000, wave_size=1000, canary_n=20, job_size=300,
                                catalog="gs://b/cat.json", ref_fa="gs://b/ref.fa", ref_fai="gs://b/ref.fa.fai")
        registry.write_run_json(paths, {"kind": "campaign", "campaign": cfg.to_dict()})
        campaign.write_samples(paths, _records(2000))
        makers = lambda c: eh.campaign_builders(c, project="proj", sa="sa@x")
        run = lambda: campaign.tick(paths, project="proj", code_version="cv1", pipeline_short="eh",
                                    make_builders=makers, log=lambda m: None)

        plan1, sub1, _, _ = run()
        assert plan1.kind == "canary" and len(sub1) == 1 and batch[sub1[0]]["n"] == 20
        assert campaign.lock_uri(paths) not in store, "lock released"
        plan2, sub2, _, _ = run()
        assert plan2.action == "wait" and not sub2

        batch[sub1[0]].update(state="SUCCEEDED", fail={3})  # 1/20 = 5%, not above the limit
        plan3, sub3, state, _ = run()
        assert plan3.kind == "wave" and len(sub3) == 4, (plan3.explain(), sub3)  # 1000 samples in 300-task jobs
        canary = state["jobs"][sub1[0]]
        assert canary["cost"]["n_vms"] == 10 and canary["cost"]["vm_usd"] > 0
        assert campaign.task_detail_uri(paths, sub1[0]) in store

        # Lose the jobs.json entry for one job: sync recovers it from Batch labels + TASKS_TSV.
        st = json.loads(store[campaign.jobs_uri(paths)])
        lost_id = sub3[0]
        lost_samples = st["jobs"].pop(lost_id)["sample_ids"]
        store[campaign.jobs_uri(paths)] = json.dumps(st)
        for j in sub3:
            batch[j].update(state="SUCCEEDED", fail=set(range(0, 300, 10)) if j == sub3[1] else set())
        plan4, sub4, state, _ = run()
        assert state["jobs"][lost_id]["sample_ids"] == lost_samples
        assert plan4.kind == "wave" and len(plan4.sample_ids) == 980, plan4.explain()

        # Stale lock from a dead notebook is broken; a fresh one blocks.
        store[campaign.lock_uri(paths)] = json.dumps({"holder": "old", "ts": 0})
        run()
        store[campaign.lock_uri(paths)] = json.dumps({"holder": "other", "ts": __import__("time").time()})
        try:
            run()
            raise AssertionError("expected LockHeld")
        except campaign.LockHeld:
            pass
    finally:
        for name, fn in saved.items():
            setattr(gcs_mod, name, fn)
        submit.list_jobs, submit.list_tasks, submit.submit_job = saved_submit
    print("campaign tick ok")


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
    test_fused_job_json()
    test_campaign_plan()
    test_job_cost_model()
    test_campaign_tick_end_to_end()
    print("all ok")


if __name__ == "__main__":
    main()
