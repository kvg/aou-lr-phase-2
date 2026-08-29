#!/usr/bin/env python3
"""Restore setup-cell tails removed by patch_notebook_bootstrap.py."""

from __future__ import annotations

import json
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]
NOTEBOOKS = REPO / "notebooks" / "terra"

BOOTSTRAP = '''from pathlib import Path
import os
import subprocess
import sys

# Bootstrap scripts/ from $WORKSPACE_BUCKET/scripts/ when not on the VM.
for _d in (Path.cwd() / "scripts", Path.cwd().parent / "scripts", Path.cwd().parent.parent / "scripts"):
    if (_d / "terra_notebook.py").is_file():
        sys.path.insert(0, str(_d.resolve()))
        break
else:
    _bucket = os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")
    if not _bucket:
        raise FileNotFoundError(
            "scripts/ not found locally and WORKSPACE_BUCKET is unset. "
            "Upload scripts/ to gs://WORKSPACE/scripts/."
        )
    _dest = (Path.cwd() / "scripts").resolve()
    _dest.mkdir(parents=True, exist_ok=True)
    subprocess.check_call(
        ["gsutil", "-m", "rsync", "-r", f"{_bucket}/scripts/", str(_dest) + "/"]
    )
    sys.path.insert(0, str(_dest))

'''


def set_code_cell(nb_path: Path, cell_idx: int, source: str) -> None:
    nb = json.loads(nb_path.read_text())
    nb["cells"][cell_idx]["source"] = source
    nb["cells"][cell_idx]["outputs"] = []
    nb["cells"][cell_idx]["execution_count"] = None
    nb_path.write_text(json.dumps(nb, indent=1) + "\n")
    print(f"updated {nb_path.name} cell {cell_idx}")


TAILS: dict[str, tuple[int, str]] = {
    "tractor_00_cov_rebuild_source.ipynb": (
        1,
        BOOTSTRAP
        + '''from __future__ import annotations

from terra_notebook import init_notebook

SCRIPTS = init_notebook("workspace_paths.py")
from workspace_paths import data_root

import json
import re
import tarfile

import numpy as np
import pandas as pd

ROOT = data_root()

RESOURCES = ROOT / "resources"
LEGACY_RESOURCES = RESOURCES / "legacy_covariates"

ONT_CSV = RESOURCES / "Multiomics_person-level_sample_counts_CDRv9 - lrWGS.phase2.ONT_v9.csv"
PACBIO_CSV = RESOURCES / "Multiomics_person-level_sample_counts_CDRv9 - lrWGS.phase2_PacBio_v9.csv"
FINAL_CSV = RESOURCES / (
    "Multiomics_person-level_sample_counts_CDRv9 - "
    "lrWGS.phase2.sample.final_set.2_Steve w v9 flag.csv"
)
MULTI_CSV = RESOURCES / "Multiomics_person-level_sample_counts_CDRv9 - multiomics_lists_CDRv9.csv"
PART_CSV = RESOURCES / "part.csv.gz"
EXTRACTION_TSV = RESOURCES / "lr_dna_extraction_method_cdrv9.tsv"
TRGT_TABLE = RESOURCES / "trgt_table.txt"
PCA_ARCHIVE = RESOURCES / "pc.log.gz"
FULL_PCA_TSV = RESOURCES / "full.tsv.gz"
CDR_MEMBERSHIP_TSV = RESOURCES / "lr.tsv"
PHENOTYPE_CSV = RESOURCES / "AoU_Phase2_Phenotype.csv.gz"
SV_COUNTS_TAR = RESOURCES / "counts.tar.gz"
SV_COLUMNS_STRUCTURE = RESOURCES / "columns_structure.txt"
ANCESTRY_CSV = LEGACY_RESOURCES / "anc_df.csv.gz"
PEDIGREE_CSV = LEGACY_RESOURCES / "aou_phase2.ped"
INTEGRATEDCALL_TSV = LEGACY_RESOURCES / "Integratedcall-GRCh38.tsv"
SV_SENS07_JSON = LEGACY_RESOURCES / "sv_sens_07_counts.json"
HA_RESEARCH_IDS_CSV = LEGACY_RESOURCES / "ha_research_ids.csv"
LR_1027_SEX_CSV = LEGACY_RESOURCES / "lr_1027_sex_at_birth.csv"
HA_PARTIAL_SAMPLE_CSV = ROOT / "ha_partial_sample.csv"
HA_PASSING_SAMPLE_CSV = LEGACY_RESOURCES / "ha_passing_sample.csv"
ONT_PHASE1_TSV = LEGACY_RESOURCES / "ont-sample-hg38.tsv"
MERGED_ALL_CSV = LEGACY_RESOURCES / "merged_all_df.csv.gz"

OUT_CSV = ROOT / "covariates.source_rebuilt.csv.gz"
OUT_DICT = ROOT / "covariates.source_rebuilt.data_dictionary.tsv"
V3_CSV = RESOURCES / "covariates.v3.csv.gz"
CONSTRUCTION_USED_V3 = False

TECH_COLS = [
    "participant",
    "long_read (1=final releasable in v9)",
    "biobank_id",
    "sex",
    "coverage",
    "platform",
    "PacBioMethylationCaller",
    "GC",
    "technology",
    "is_AIAN",
    "withdrew",
]
MULTI_COLS = [
    "research_id",
    "long_read (1=meet QC requirement)",
    "long_read phase",
    "long_read (1=final releasable in v9)",
    "RNASeq (1=meet QC requirement)",
    "RNASeq(1=final releasable)",
    "Proteomics (1=meet QC requirement)",
    "Proteomics(1=final releasable)",
    "Exposomics (1=releasable in v9)",
    "rid_no_srWGS",
]
PART_COLS = ["person_id", "age_at_cdr", "zip3_as_string", "sex_at_birth", "has_ehr_data"]
PC_COLS = [f"PC{i}" for i in range(1, 33)]
EXTRACTION_COLS = ["research_id", "extraction_method"]

print("ROOT:", ROOT)
print("RESOURCES:", RESOURCES)
for p in [
    ONT_CSV,
    PACBIO_CSV,
    FINAL_CSV,
    MULTI_CSV,
    PART_CSV,
    EXTRACTION_TSV,
    TRGT_TABLE,
    PCA_ARCHIVE,
    FULL_PCA_TSV,
    CDR_MEMBERSHIP_TSV,
    PHENOTYPE_CSV,
    SV_COUNTS_TAR,
    SV_COLUMNS_STRUCTURE,
    ANCESTRY_CSV,
    PEDIGREE_CSV,
    INTEGRATEDCALL_TSV,
    SV_SENS07_JSON,
    HA_RESEARCH_IDS_CSV,
    LR_1027_SEX_CSV,
    HA_PARTIAL_SAMPLE_CSV,
    HA_PASSING_SAMPLE_CSV,
    ONT_PHASE1_TSV,
    MERGED_ALL_CSV,
]:
    assert p.exists(), p
    print("  ok", p.name)
''',
    ),
    "tractor_02_qc_results.ipynb": (
        1,
        BOOTSTRAP
        + '''from terra_notebook import init_notebook

SCRIPTS = init_notebook(
    "workspace_paths.py",
    "compare_calibration.py",
    "plot_tractor_results.py",
)
from workspace_paths import data_root

ROOT = data_root()
COMPARE_PY = SCRIPTS / "compare_calibration.py"
PLOT_PY = SCRIPTS / "plot_tractor_results.py"

# Result directories (local). Leave unset / empty to skip a model.
TRACTOR_LIMITED_DIR = Path(os.environ.get("TRACTOR_LIMITED_DIR", "results/tractor_limited"))
TRACTOR_FULL_DIR = Path(os.environ.get("TRACTOR_FULL_DIR", "results/tractor_full"))
SAIGE_LIMITED_DIR = Path(os.environ.get("SAIGE_LIMITED_DIR", "results/saige_limited"))
SAIGE_FULL_DIR = Path(os.environ.get("SAIGE_FULL_DIR", "results/saige_full"))
NULL_META_DIR = Path(os.environ.get("SAIGE_NULL_META_DIR", "results/saige_null_meta"))

OUT_DIR = Path(os.environ.get("CALIBRATION_QC_DIR", "calibration_qc"))
OUT_DIR.mkdir(parents=True, exist_ok=True)

ws = os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")


def sh(cmd: str) -> None:
    print(cmd)
    rc = get_ipython().system(cmd)
    if rc:
        raise RuntimeError(f"command failed with exit code {rc}: {cmd}")


print("OUT_DIR:", OUT_DIR.resolve())
for label, d in [
    ("tractor_limited", TRACTOR_LIMITED_DIR),
    ("tractor_full", TRACTOR_FULL_DIR),
    ("saige_limited", SAIGE_LIMITED_DIR),
    ("saige_full", SAIGE_FULL_DIR),
]:
    n = len(list(d.glob("*.tsv"))) if d.exists() else 0
    print(f"  {label}: {d} ({n} tsv)")
''',
    ),
    "tractor_03_cov_summarize.ipynb": (
        1,
        BOOTSTRAP
        + '''from __future__ import annotations

from terra_notebook import init_notebook

SCRIPTS = init_notebook("workspace_paths.py")
from workspace_paths import data_root

import warnings

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

warnings.filterwarnings("ignore", category=FutureWarning)
pd.set_option("display.max_columns", 80)
pd.set_option("display.max_rows", 80)
pd.set_option("display.width", 140)
sns.set_theme(style="whitegrid", context="notebook")
plt.rcParams.update({
    "figure.dpi": 120,
    "savefig.dpi": 200,
    "font.size": 11,
    "axes.titlesize": 13,
    "axes.labelsize": 11,
    "legend.fontsize": 9,
})

ROOT = data_root()
COV_CSV = ROOT / "covariates.source_rebuilt.csv.gz"
COV_DICT = ROOT / "covariates.source_rebuilt.data_dictionary.tsv"
PHENOTYPE_CSV = ROOT / "resources" / "AoU_Phase2_Phenotype.csv.gz"
OUT_DIR = ROOT / "summaries"
FIG_DIR = OUT_DIR / "figures"
TAB_DIR = OUT_DIR / "tables"
MS_DIR = OUT_DIR / "manuscript"

for d in [OUT_DIR, FIG_DIR, TAB_DIR, MS_DIR]:
    d.mkdir(parents=True, exist_ok=True)

MIN_CELL = 20
ARTIFACTS: list[Path] = []

print("ROOT:", ROOT)
print("OUT_DIR:", OUT_DIR)
''',
    ),
    "tractor_04_table1_cohort_summary.ipynb": (
        1,
        BOOTSTRAP
        + '''from __future__ import annotations

from terra_notebook import init_notebook

SCRIPTS = init_notebook("workspace_paths.py")
from workspace_paths import data_root

import numpy as np
import pandas as pd

ROOT = data_root()
RESOURCES = ROOT / "resources"
LEGACY = RESOURCES / "legacy_covariates"

COV_CSV = ROOT / "covariates.source_rebuilt.csv.gz"
INTEGRATEDCALL_TSV = LEGACY / "Integratedcall-GRCh38.tsv"
MERGED_ALL_CSV = LEGACY / "merged_all_df.csv.gz"
ONT_PHASE1_TSV = LEGACY / "ont-sample-hg38.tsv"
PEDIGREE_CSV = LEGACY / "aou_phase2.ped"

OUT_DIR = ROOT / "summaries" / "manuscript"
OUT_DIR.mkdir(parents=True, exist_ok=True)
OUT_TSV = OUT_DIR / "table1_cohort_summary.tsv"
OUT_MD = OUT_DIR / "table1_cohort_summary.md"

HIGH_PASS_N = 1133

for path in [COV_CSV, INTEGRATEDCALL_TSV, MERGED_ALL_CSV, ONT_PHASE1_TSV, PEDIGREE_CSV]:
    assert path.exists(), path

print("ROOT:", ROOT)
print("OUT_DIR:", OUT_DIR)
''',
    ),
    "tractor_05_pca_deepvariant_long_read.ipynb": (
        1,
        BOOTSTRAP
        + '''from __future__ import annotations

from terra_notebook import init_notebook

SCRIPTS = init_notebook("workspace_paths.py")
from workspace_paths import data_root

import json
import os

import pandas as pd

try:
    display
except NameError:
    def display(value):
        print(value)

ROOT = data_root()
WORKSPACE_BUCKET = os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")
RUN_LABEL = os.environ.get("PCA_RUN_LABEL", "deepvariant_lr_v2")

if WORKSPACE_BUCKET:
    _bucket = (
        WORKSPACE_BUCKET
        if WORKSPACE_BUCKET.startswith("gs://")
        else f"gs://{WORKSPACE_BUCKET}"
    )
    COV_URI = os.environ.get(
        "PCA_COV_URI", f"{_bucket}/covariates/covariates.source_rebuilt.csv.gz"
    )
    CHROM_MANIFEST_URI = os.environ.get(
        "PCA_CHROM_MANIFEST_URI", f"{_bucket}/manifests/GL_INTERVAL_set.tsv"
    )
    OUT_DIR = os.environ.get("PCA_OUT_DIR", f"{_bucket}/pca/{RUN_LABEL}")
    EXISTING_MT = os.environ.get("PCA_EXISTING_MT", f"{_bucket}/mt/deepvariant_joint.mt")
    EXISTING_VDS = os.environ.get("PCA_EXISTING_VDS", f"{_bucket}/vds/deepvariant_joint.vds")
else:
    COV_URI = os.environ.get("PCA_COV_URI", str(ROOT / "covariates.source_rebuilt.csv.gz"))
    CHROM_MANIFEST_URI = os.environ.get(
        "PCA_CHROM_MANIFEST_URI", str(ROOT / "manifests" / "GL_INTERVAL_set.tsv")
    )
    OUT_DIR = os.environ.get("PCA_OUT_DIR", str(ROOT / "pca" / RUN_LABEL))
    EXISTING_MT = os.environ.get("PCA_EXISTING_MT", "")
    EXISTING_VDS = os.environ.get("PCA_EXISTING_VDS", "")

INPUT_MODE = os.environ.get("PCA_INPUT_MODE", "vcf_list")
CHROM_ID_COLUMN = "entity:GL_INTERVAL_set_id"
CHROM_URI_COLUMN = "VCF"
CHROM_IDX_COLUMN = "VCF_idx"
VCF_URIS: list[str] | None = None
AUTOSOMES_ONLY = True

N_PCS = int(os.environ.get("PCA_N_PCS", "32"))
MIN_AF = 0.01
MAX_AF = 0.99
MIN_VARIANT_CALL_RATE = 0.98
MIN_SAMPLE_CALL_RATE = 0.98
LD_R2 = 0.1
LD_BP_WINDOW = 500_000
LD_MEMORY_PER_CORE = os.environ.get("PCA_LD_MEMORY_PER_CORE", "1g")
PASS_ONLY = True
MIN_POPULATION_N = 100
POPULATION_LABEL = "ancestry_pred_other"
ALLOW_PARTIAL_VCF_OVERLAP = True
RUN_PIPELINE = os.environ.get("PCA_RUN_PIPELINE", "").lower() in {"1", "true", "yes"}

CHECKPOINT_MT = f"{OUT_DIR}/checkpoints/qc_for_pca.mt"
GLOBAL_PCS_TSV = f"{OUT_DIR}/global_pcs.tsv"
POP_PCS_TSV = f"{OUT_DIR}/population_pcs.tsv"
METADATA_JSON = f"{OUT_DIR}/run_metadata.json"

print("ROOT:", ROOT)
print("WORKSPACE_BUCKET:", WORKSPACE_BUCKET or "(local)")
print("OUT_DIR:", OUT_DIR)
print("INPUT_MODE:", INPUT_MODE)
print("RUN_PIPELINE:", RUN_PIPELINE)
''',
    ),
    "tractor_08_qc_grm.ipynb": (
        1,
        BOOTSTRAP
        + '''from terra_notebook import init_notebook

SCRIPTS = init_notebook("workspace_paths.py", "plot_grm_qc.R")
from workspace_paths import data_root

import shutil
import subprocess

ROOT = data_root()
ws = os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")

GRM_SPARSE_RDS = os.environ.get("TRACTOR_GRM_SPARSE_RDS", "")
GRM_BANDS = os.environ.get("TRACTOR_GRM_BANDS", "")
GRM_HIST = os.environ.get("TRACTOR_GRM_HIST", "")
GRM_CLOSE = os.environ.get("TRACTOR_GRM_CLOSE", "")
COV_GCS = os.environ.get(
    "TRACTOR_COVARIATES_GCS",
    f"{ws}/covariates/covariates.source_rebuilt.csv.gz" if ws else "",
)

OUT = Path(os.environ.get("GRM_QC_OUT", "grm_qc"))
OUT.mkdir(parents=True, exist_ok=True)
LOCAL = Path(os.environ.get("GRM_QC_WORK", "grm_qc_inputs"))
LOCAL.mkdir(parents=True, exist_ok=True)


def pull(uri: str, name: str) -> Path | None:
    if not uri:
        return None
    dest = LOCAL / name
    if dest.exists():
        return dest
    if uri.startswith("gs://"):
        subprocess.check_call(["gsutil", "cp", uri, str(dest)])
    else:
        shutil.copy2(uri, dest)
    return dest


assert GRM_SPARSE_RDS, "Set TRACTOR_GRM_SPARSE_RDS to the MakeGRM grm_sparse.rds URI"
rds = pull(GRM_SPARSE_RDS, "grm_sparse.rds")
bands = pull(GRM_BANDS, "grm_relationship_bands.tsv")
hist = pull(GRM_HIST, "grm_kinship_histogram.tsv")
close = pull(GRM_CLOSE, "grm_close_pairs.tsv")
cov = pull(COV_GCS, "covariates.source_rebuilt.csv.gz") if COV_GCS else None
print("inputs:", {"rds": rds, "bands": bands, "hist": hist, "close": close, "cov": cov})
''',
    ),
    "sv_02_stage_sample_ancestry.ipynb": (
        1,
        BOOTSTRAP
        + '''from terra_notebook import init_notebook

SCRIPTS = init_notebook("workspace_paths.py", "sv_site_utils.py")
from workspace_paths import data_root

import gzip
import shutil
import subprocess

import pandas as pd
from sv_site_utils import order_samples_by_ancestry  # noqa: E402

ROOT = data_root()

WORK = Path(os.environ.get("SAMPLE_ANCESTRY_WORK", Path.cwd() / "sample_ancestry_work")).resolve()
WORK.mkdir(parents=True, exist_ok=True)

WORKSPACE_BUCKET = os.environ.get("WORKSPACE_BUCKET", "").rstrip("/")
GCS_PREFIX = os.environ.get(
    "SAMPLE_ANCESTRY_GCS_PREFIX",
    f"{WORKSPACE_BUCKET}/metadata" if WORKSPACE_BUCKET else "",
)
COV_GCS = os.environ.get(
    "TRACTOR_COVARIATES_GCS",
    f"{WORKSPACE_BUCKET}/covariates/covariates.source_rebuilt.csv.gz"
    if WORKSPACE_BUCKET
    else "",
)
COV_CSV = Path(os.environ.get("SAMPLE_ANCESTRY_COV_CSV", ROOT / "covariates.source_rebuilt.csv.gz"))

PHASE1_MAIN_VCF = os.environ.get("SV_PHASE1_MAIN_VCF", "").strip()
PHASE2_MAIN_VCF = os.environ.get("SV_PHASE2_MAIN_VCF", "").strip()

CANONICAL = {"eas", "amr", "eur", "sas", "afr", "oth"}
TO_OTH = {"", "nan", "none", "unknown", "other", "mid", "oth"}

print("ROOT:", ROOT)
print("COV_CSV:", COV_CSV)
print("COV_GCS:", COV_GCS or "(unset)")
print("PHASE1_MAIN_VCF:", PHASE1_MAIN_VCF or "(lr_phase filter)")
print("PHASE2_MAIN_VCF:", PHASE2_MAIN_VCF or "(lr_phase filter)")
print("GCS_PREFIX:", GCS_PREFIX or "(set WORKSPACE_BUCKET or SAMPLE_ANCESTRY_GCS_PREFIX)")
''',
    ),
    "sv_03_manuscript_stats.ipynb": (
        1,
        BOOTSTRAP
        + '''from terra_notebook import init_notebook

SCRIPTS = init_notebook(
    "workspace_paths.py",
    "merge_manuscript_counts.py",
    "plot_discovery.py",
)
from workspace_paths import sv_output_dir

import pandas as pd

OUT = sv_output_dir()

PHASE1_COUNTS = OUT / "aou_lr_phase1.manuscript_counts.tsv"
PHASE2_COUNTS = OUT / "aou_lr_phase2.manuscript_counts.tsv"
PHASE2_DISCOVERY = OUT / "aou_lr_phase2.discovery.tsv"
PHASE2_DISCOVERY_REGION = OUT / "aou_lr_phase2.discovery.region.tsv"
PHASE2_DISCOVERY_CADD = OUT / "aou_lr_phase2.discovery.cadd.tsv"
OUTDIR = OUT / "figures"
OUTDIR.mkdir(parents=True, exist_ok=True)

print("OUT:", OUT)
print("OUTDIR:", OUTDIR)
''',
    ),
}


def fix_tractor_cov_explore() -> None:
    path = NOTEBOOKS / "tractor_cov_explore.ipynb"
    nb = json.loads(path.read_text())
    nb["cells"][0]["cell_type"] = "markdown"
    nb["cells"][0]["source"] = ["# Ad hoc covariate exploration\n"]
    load_cell = (
        BOOTSTRAP
        + '''from __future__ import annotations

from terra_notebook import init_notebook

SCRIPTS = init_notebook("workspace_paths.py")
from workspace_paths import data_root

import sys

import matplotlib.pyplot as plt
import pandas as pd
import seaborn as sns

pd.set_option("display.max_columns", 80)
pd.set_option("display.max_rows", 80)
pd.set_option("display.width", 140)
sns.set_theme(style="whitegrid", context="notebook")

ROOT = data_root()
COV_CSV = ROOT / "covariates.source_rebuilt.csv.gz"
COV_DICT = ROOT / "covariates.source_rebuilt.data_dictionary.tsv"
print("ROOT:", ROOT)
print("CSV:", COV_CSV)
'''
    )
    nb["cells"][2]["source"] = load_cell
    nb["cells"][2]["outputs"] = []
    nb["cells"][2]["execution_count"] = None
    path.write_text(json.dumps(nb, indent=1) + "\n")
    print("updated tractor_cov_explore.ipynb")


def main() -> None:
    for name, (idx, tail) in TAILS.items():
        path = NOTEBOOKS / name
        if not path.is_file():
            print(f"skip missing {name}")
            continue
        set_code_cell(path, idx, tail)
    fix_tractor_cov_explore()


if __name__ == "__main__":
    main()
