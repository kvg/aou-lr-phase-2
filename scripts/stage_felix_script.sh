#!/usr/bin/env bash
# Stage ONE file from felix/scripts/ to the workspace bucket, verifying both ends.
#
# Use this instead of stage_felix_scripts.sh during a run. The wholesale stager
# uploads eight files, and fit_felix_null.R and run_felix_step2.R are
# content-hashed File inputs to FitFelixNull and RunFelixStep2 -- re-staging a
# drifted copy of either invalidates every cached null and every step-2 shard.
#
# It also checks afterwards. A gsutil cp run from a stale checkout reports
# success and changes nothing, which is how fd6f519 of summarize_felix_results.py
# stayed in the bucket through an upload that looked like it worked.
#
# Usage:  ./scripts/stage_felix_script.sh summarize_felix_results.py
set -euo pipefail

NAME="${1:?usage: stage_felix_script.sh <filename under felix/scripts/>}"
: "${WORKSPACE_BUCKET:?set WORKSPACE_BUCKET}"
REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SRC="${REPO_ROOT}/felix/scripts/${NAME}"
DEST="${WORKSPACE_BUCKET%/}/felix/scripts/${NAME}"

[[ -f "$SRC" ]] || { echo "no such file: $SRC" >&2; exit 1; }

md5of() { md5sum "$1" 2>/dev/null | cut -d' ' -f1 || md5 -q "$1"; }
local_md5=$(md5of "$SRC")

# The checkout must be current, or you upload the same stale bytes again.
if command -v git >/dev/null && git -C "$REPO_ROOT" rev-parse --git-dir >/dev/null 2>&1; then
  git -C "$REPO_ROOT" fetch -q origin 2>/dev/null || true
  behind=$(git -C "$REPO_ROOT" rev-list --count HEAD..origin/main 2>/dev/null || echo 0)
  dirty=$(git -C "$REPO_ROOT" status --porcelain -- "felix/scripts/${NAME}")
  echo "checkout: HEAD $(git -C "$REPO_ROOT" rev-parse --short HEAD), ${behind} commit(s) behind origin/main"
  [[ -n "$dirty" ]] && echo "  NOTE: ${NAME} is modified relative to HEAD"
  if [[ "$behind" != "0" ]]; then
    echo "  refusing to stage from a checkout behind origin/main; run 00_sync_repo or git pull" >&2
    exit 1
  fi
fi

echo "local  $local_md5  $SRC"
remote_before=$(gsutil cat "$DEST" 2>/dev/null | md5sum | cut -d' ' -f1 || echo "<absent>")
echo "remote $remote_before  $DEST"
if [[ "$remote_before" == "$local_md5" ]]; then
  echo "already current, nothing to do"; exit 0
fi

gsutil cp "$SRC" "$DEST"

remote_after=$(gsutil cat "$DEST" | md5sum | cut -d' ' -f1)
echo "remote $remote_after  (after upload)"
if [[ "$remote_after" != "$local_md5" ]]; then
  echo "UPLOAD DID NOT TAKE: remote still $remote_after, expected $local_md5" >&2
  exit 1
fi
echo "staged ${NAME}; the tasks consuming it will cache-miss and re-run, nothing else will"
