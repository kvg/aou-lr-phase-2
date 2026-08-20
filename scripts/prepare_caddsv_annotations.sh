#!/usr/bin/env bash
# Download and archive CADD-SV v2 annotations once. Upload the resulting tar.gz
# to the workspace/reference bucket and provide it as caddsv_annotations_tar.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
IMAGE="${SV_ANNOTATION_IMAGE:-us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.5}"
OUT="${1:-$ROOT/caddsv-v2.0-annotations.tar.gz}"

mkdir -p "$(dirname "$OUT")"
available_kb="$(df -Pk "$(dirname "$OUT")" | awk 'NR==2 {print $4}')"
if [[ "$available_kb" -lt 60000000 ]]; then
  echo "Need at least 60 GB free where the Docker engine stores data to package CADD-SV annotations." >&2
  echo "Run this on a larger VM or use a Docker engine with a larger disk, then upload the archive to GCS." >&2
  exit 1
fi

TMP="$(mktemp -d)"
trap 'rm -rf "$TMP"' EXIT

docker run --rm --platform linux/amd64 \
  -v "$TMP:/work" \
  "$IMAGE" \
  bash -c 'caddsv get annotations --annotations-dir /work/annotations && tar -C /work -czf /work/caddsv-annotations.tar.gz annotations'

mv "$TMP/caddsv-annotations.tar.gz" "$OUT"
echo "Wrote $OUT"
