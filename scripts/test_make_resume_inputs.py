"""felix/scripts/make_resume_inputs.py against a stub gsutil serving a fake run listing."""
import json
import stat
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).resolve().parents[1] / "felix" / "scripts" / "make_resume_inputs.py"
ROOT = "gs://bkt/submissions/sub1/FelixPilot/wf1"
PHENOS = ["CV_401", "GI_520", "NS_300"]

STUB = """#!{python}
import json, os, sys
fx = json.load(open(os.environ["STUB_FIXTURE"]))
if sys.argv[1] == "ls":
    assert sys.argv[2].endswith("/**"), sys.argv
    pre = sys.argv[2][:-3]  # strip the trailing /**
    hits = [o for o in fx["objects"] if o.startswith(pre + "/")]
    if not hits:
        sys.stderr.write("CommandException: One or more URLs matched no objects.\\n"); sys.exit(1)
    print("\\n".join(hits))
elif sys.argv[1] == "cat":
    print(fx["cat"][sys.argv[2]], end="")
"""


def objects(extra=(), skip=()):
    objs = [f"{ROOT}/call-Pack/felixla_packed.tar.gz", f"{ROOT}/call-Pack/stderr",
            f"{ROOT}/call-MakeGRM/saige_sparseGRM.mtx", f"{ROOT}/call-MakeGRM/saige_sparseGRM.sampleIDs.txt",
            f"{ROOT}/call-MakeGRM/saige_plink.bed"]
    for i, p in enumerate(PHENOS):
        for suffix in ("null.rda", "varianceRatio.txt", "samples.txt", "null_meta.tsv", "pheno.tsv"):
            objs.append(f"{ROOT}/call-Null/shard-{i}/{p}.{suffix}")
        objs += [f"{ROOT}/call-Null/shard-{i}/rc", f"{ROOT}/call-Null/shard-{i}/stderr"]
    objs += list(extra)
    return [o for o in objs if o not in set(skip)]


@pytest.fixture
def env(tmp_path):
    stub = tmp_path / "gsutil"
    stub.write_text(STUB.replace("{python}", sys.executable))
    stub.chmod(stub.stat().st_mode | stat.S_IEXEC)
    base = tmp_path / "base.json"
    base.write_text(json.dumps({
        "_comment": "x", "FelixPilot.pheno_cov": "gs://bkt/pheno_cov.tsv",
        "FelixPilot.run_felix_step2_script": "gs://bkt/felix/scripts/run_felix_step2.R",
        "FelixPilot.summarize_script": "gs://bkt/felix/scripts/summarize_felix_results.py",
        "FelixPilot.chrom": "chr22", "FelixPilot.num_ancs": 5, "FelixPilot.min_mac": 50,
        "FelixPilot.docker": "img:0.2.0", "FelixPilot.fit_felix_null_script": "not copied"}))
    return tmp_path, stub, base


def run(env, objs, limit=0, phenos_arg=None):
    tmp, stub, base = env
    fx = tmp / "fx.json"
    fx.write_text(json.dumps({"objects": objs, "cat": {"gs://bkt/sel.txt": "\n".join(PHENOS) + "\n"}}))
    out = tmp / "out.json"
    cmd = [sys.executable, str(SCRIPT), "--run-root", ROOT + "/", "--phenotypes-file", phenos_arg or "gs://bkt/sel.txt",
           "--base-config", str(base), "--out", str(out), "--gsutil", str(stub), "--limit", str(limit)]
    proc = subprocess.run(cmd, capture_output=True, text=True, env={"STUB_FIXTURE": str(fx), "PATH": ""})
    return proc, out


def test_all_phenotypes_in_order(env):
    proc, out = run(env, objects())
    assert proc.returncode == 0, proc.stderr
    d = json.loads(out.read_text())
    p = "FelixPilotStep2Resume."
    assert d[p + "phenotypes"] == PHENOS
    assert d[p + "null_rdas"] == [f"{ROOT}/call-Null/shard-{i}/{x}.null.rda" for i, x in enumerate(PHENOS)]
    assert d[p + "variance_ratios"][1].endswith("shard-1/GI_520.varianceRatio.txt")
    assert d[p + "samples_used"][2].endswith("shard-2/NS_300.samples.txt")
    assert d[p + "packed_tar"] == f"{ROOT}/call-Pack/felixla_packed.tar.gz"
    assert d[p + "sparse_grm_mtx"].endswith("saige_sparseGRM.mtx")
    assert d[p + "sparse_grm_sample_ids"].endswith("saige_sparseGRM.sampleIDs.txt")
    assert d[p + "pheno_cov"] == "gs://bkt/pheno_cov.tsv" and d[p + "num_ancs"] == 5
    assert p + "fit_felix_null_script" not in d and "_comment" not in str(d.keys())


def test_limit_and_local_phenotype_file(env):
    sel = env[0] / "sel.txt"
    sel.write_text("\n".join(PHENOS) + "\n")
    proc, out = run(env, objects(), limit=2, phenos_arg=str(sel))
    assert proc.returncode == 0, proc.stderr
    d = json.loads(out.read_text())
    assert d["FelixPilotStep2Resume.phenotypes"] == PHENOS[:2]
    assert len(d["FelixPilotStep2Resume.null_rdas"]) == 2


def test_duplicate_attempt_is_refused(env):
    dup = f"{ROOT}/call-Null/shard-1/attempt-2/GI_520.null.rda"
    proc, _ = run(env, objects(extra=[dup]))
    assert proc.returncode != 0
    assert "shard-1" in proc.stderr and "found 2" in proc.stderr and dup in proc.stderr


def test_missing_shard_and_missing_file_are_refused(env):
    proc, _ = run(env, objects(skip=[f"{ROOT}/call-Null/shard-2/NS_300.samples.txt"]))
    assert proc.returncode != 0 and "NS_300.samples.txt" in proc.stderr and "found 0" in proc.stderr
    gone = [o for o in objects() if "/shard-2/" not in o]
    proc, _ = run(env, gone)
    assert proc.returncode != 0 and "shard-2" in proc.stderr
