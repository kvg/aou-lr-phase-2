#!/usr/bin/env bash
# Summarize Cromwell/Terra task timing from a downloaded metadata JSON directory.
#
# Usage:
#   ./scripts/summarize_terra_timing.sh /path/to/metadata/
#   ./scripts/summarize_terra_timing.sh gs://bucket/submission/metadata/
#
# Pull metadata first (Terra Job Manager → export, or):
#   gsutil -m rsync -r gs://.../submissions/<id>/ metadata_old/
#
# Prints per-call-name wall times and a Score-vs-total breakdown when present.

set -euo pipefail

SRC="${1:?usage: summarize_terra_timing.sh METADATA_DIR_OR_GCS}"

WORKDIR="$(mktemp -d)"
trap 'rm -rf "${WORKDIR}"' EXIT

if [[ "${SRC}" == gs://* ]]; then
  gsutil -m rsync -r "${SRC%/}/" "${WORKDIR}/"
  META_ROOT="${WORKDIR}"
else
  META_ROOT="${SRC}"
fi

python3 - <<PY
import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

root = Path("${META_ROOT}")

def parse_ts(s):
    if not s:
        return None
    s = s.replace("Z", "+00:00")
    return datetime.fromisoformat(s)

rows = []
for path in sorted(root.rglob("*.json")):
    try:
        meta = json.loads(path.read_text())
    except Exception:
        continue
    if not isinstance(meta, dict):
        continue
    name = meta.get("callRoot") or meta.get("name") or path.stem
    start = parse_ts(meta.get("startTime") or meta.get("start"))
    end = parse_ts(meta.get("endTime") or meta.get("end"))
    if not start or not end:
        continue
    sec = (end - start).total_seconds()
    runtime = meta.get("runtimeAttributes") or {}
    cpu = runtime.get("cpu") or runtime.get("cpus") or "?"
    mem = runtime.get("memory") or runtime.get("memoryInGB") or "?"
    rows.append((name, sec, cpu, mem, str(path.relative_to(root))))

if not rows:
    print("No timing JSON found under", root)
    raise SystemExit(1)

by_call = defaultdict(list)
for name, sec, cpu, mem, rel in rows:
    # collapse scatter indices: TractorMixPilot.Score-7 -> TractorMixPilot.Score
    base = name.split("-")[0] if ".Score-" in name or ".FitNull-" in name else name
    by_call[base].append((sec, cpu, mem, name))

print(f"{'call':<40} {'n':>4} {'total_h':>8} {'mean_min':>9} {'max_min':>9} {'cpu':>4} {'mem':>8}")
print("-" * 85)
for call in sorted(by_call):
    items = by_call[call]
    total = sum(x[0] for x in items)
    mean = total / len(items)
    mx = max(x[0] for x in items)
    cpu = items[0][1]
    mem = items[0][2]
    print(f"{call:<40} {len(items):4d} {total/3600:8.2f} {mean/60:9.1f} {mx/60:9.1f} {str(cpu):>4} {str(mem):>8}")

score = [sec for name, sec, _, _, _ in rows if ".Score" in name]
fit = [sec for name, sec, _, _, _ in rows if ".FitNull" in name]
if score:
    print()
    print(f"Score tasks: n={len(score)} total={sum(score)/3600:.2f}h mean={sum(score)/len(score)/60:.1f}min")
if fit:
    print(f"FitNull tasks: n={len(fit)} total={sum(fit)/3600:.2f}h mean={sum(fit)/len(fit)/60:.1f}min")
if score and fit:
    print(f"Score / FitNull wall-time ratio: {sum(score)/sum(fit):.2f}x")
PY
