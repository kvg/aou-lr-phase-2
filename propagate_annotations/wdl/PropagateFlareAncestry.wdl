version 1.0

# Copy FLARE FORMAT/AN1,AN2 onto interstitial sites of a target (phased) VCF.
#
# FLARE markers are a sparse LAI subset. Every target record keeps its GT and
# inherits AN1/AN2 from the covering FLARE interval (FELIXla semantics).
# Both inputs are tabix/csi-indexed; the task subsets with bcftools view -r,
# left-padding FLARE to chrom:1-END so the marker before `region` is available,
# then streams the two extracts in linear time.
#
# On Terra (Pipelines API), VCF/index File inputs are localization_optional:
# Cromwell passes gs:// URIs and bcftools/htslib range-reads the region.
# See https://cromwell.readthedocs.io/en/latest/optimizations/FileLocalization/
#
#   flare_vcf            = this.model_chr_anc_vcf   (aou_lr_chrom, or FlareByPopulation anc VCF)
#   target_vcf           = phased chrom VCF (SNV/indel, SV, or joint)
#   propagate_script    = gs://$WORKSPACE_BUCKET/scripts/propagate_flare_ancestry.py
#   gcs_project          = Terra google/billing project (requester-pays)

workflow PropagateFlareAncestry {
  input {
    File flare_vcf
    File flare_vcf_index
    File target_vcf
    File target_vcf_index
    String region
    String prefix = "flare_ancestry"
    File propagate_script
    String gcs_project = ""

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-flare-bcftools:1.24"
    Float min_coverage_rate = 0.95
    Boolean fail_on_low_coverage = false
    Float min_sample_overlap = 0.0
    String format_tags = "AN1,AN2"
    Int cpu = 4
    Int memory_gb = 8
    Int preemptible = 1
    Int max_retries = 3
    Int stream_attempts = 3
    # Region extracts + output only. Full VCFs stay in GCS when streaming.
    # If Cromwell localizes (local backend), set disk_gb_multiplier ~ 2.0.
    Int disk_gb_floor = 40
    Float disk_gb_multiplier = 0.0
    Int bcftools_threads = 4
  }

  call AnnotateTarget {
    input:
      flare_vcf = flare_vcf,
      flare_vcf_index = flare_vcf_index,
      target_vcf = target_vcf,
      target_vcf_index = target_vcf_index,
      region = region,
      prefix = prefix,
      propagate_script = propagate_script,
      gcs_project = gcs_project,
      docker = docker,
      min_coverage_rate = min_coverage_rate,
      fail_on_low_coverage = fail_on_low_coverage,
      min_sample_overlap = min_sample_overlap,
      format_tags = format_tags,
      cpu = cpu,
      memory_gb = memory_gb,
      preemptible = preemptible,
      max_retries = max_retries,
      stream_attempts = stream_attempts,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier,
      bcftools_threads = bcftools_threads
  }

  output {
    File annotated_vcf = AnnotateTarget.annotated_vcf
    File annotated_vcf_index = AnnotateTarget.annotated_vcf_index
    File stats_json = AnnotateTarget.stats_json
    File stats_tsv = AnnotateTarget.stats_tsv
    File missing_tsv = AnnotateTarget.missing_tsv
  }

  meta {
    description: "Region-subset target VCF annotated with covering FLARE AN1/AN2 (linear streaming, FELIXla intervals)."
    allowNestedInputs: true
  }
}

task AnnotateTarget {
  input {
    File flare_vcf
    File flare_vcf_index
    File target_vcf
    File target_vcf_index
    String region
    String prefix
    File propagate_script
    String gcs_project
    String docker
    Float min_coverage_rate
    Boolean fail_on_low_coverage
    Float min_sample_overlap
    String format_tags
    Int cpu
    Int memory_gb
    Int preemptible
    Int max_retries
    Int stream_attempts
    Int disk_gb_floor
    Float disk_gb_multiplier
    Int bcftools_threads
  }

  Int disk_gb = ceil((size(flare_vcf, "GB") + size(target_vcf, "GB")) * disk_gb_multiplier) + disk_gb_floor

  parameter_meta {
    flare_vcf: {
      description: "FLARE ancestry VCF. Streamed from GCS on Terra PAPI; must also work as a local path.",
      localization_optional: true
    }
    flare_vcf_index: {
      description: "Tabix/CSI sibling of flare_vcf.",
      localization_optional: true
    }
    target_vcf: {
      description: "Phased target VCF (SNV/indel, SV, or joint). GT is kept; AN1/AN2 are added.",
      localization_optional: true
    }
    target_vcf_index: {
      description: "Tabix/CSI sibling of target_vcf.",
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail
    python3 -c "import pathlib,sys; pathlib.Path('region.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{region}"
    REGION="$(cat region.txt)"
    python3 - <<'PY'
from pathlib import Path
spec = Path("region.txt").read_text().strip()
if ":" not in spec:
    Path("flare_region.txt").write_text(spec + "\n")
else:
    chrom, rest = spec.split(":", 1)
    end = rest.split("-", 1)[-1]
    Path("flare_region.txt").write_text(f"{chrom}:1-{end}\n")
PY
    FLARE_REGION="$(cat flare_region.txt)"
    echo "[$(date -Is)] PropagateFlareAncestry target=${REGION} flare=${FLARE_REGION}" >&2
    echo "FLARE URI: ~{flare_vcf}" >&2
    echo "target URI: ~{target_vcf}" >&2
    df -h . || true
    bcftools --version | head -n 2 >&2 || true

    cat > gcs_auth.py <<'PY'
import json, os, subprocess, sys, urllib.request
from pathlib import Path

def fresh_token():
    for cmd in (
        ["gcloud", "auth", "print-access-token"],
        ["gcloud", "auth", "application-default", "print-access-token"],
    ):
        try:
            tok = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
            if tok:
                return tok
        except Exception:
            pass
    for url in (
        "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
        "http://169.254.169.254/computeMetadata/v1/instance/service-accounts/default/token",
    ):
        req = urllib.request.Request(url, headers={"Metadata-Flavor": "Google"})
        try:
            with urllib.request.urlopen(req, timeout=5) as resp:
                tok = json.load(resp).get("access_token") or ""
                if tok:
                    return tok
        except Exception:
            continue
    return ""

need_gcs = any(u.startswith("gs://") for u in sys.argv[1:])
token = fresh_token() if need_gcs else ""
if need_gcs and not token:
    raise SystemExit("error: bcftools needs a GCS OAuth token to read gs:// VCFs")
project = "~{gcs_project}".strip() or os.environ.get("GOOGLE_PROJECT", "")
parts = []
if project:
    parts.append(f"export GCS_REQUESTER_PAYS_PROJECT={project}")
if token:
    parts.append(f"export GCS_OAUTH_TOKEN='{token}'")
Path("gcs.env").write_text(("\n".join(parts) + "\n") if parts else "")
PY

    refresh_gcs_auth() {
      python3 gcs_auth.py "~{flare_vcf}" "~{target_vcf}"
      if [[ -s gcs.env ]]; then
        # shellcheck disable=SC1091
        source gcs.env
      fi
    }

    bcftools_view_retry() {
      local out="$1"
      shift
      local attempt=1
      local max=~{stream_attempts}
      local delay=20
      if [[ "${max}" -lt 1 ]]; then
        max=1
      fi
      while [[ "${attempt}" -le "${max}" ]]; do
        refresh_gcs_auth
        echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) bcftools view attempt ${attempt}/${max} -> ${out}" >&2
        rm -f "${out}"
        if bcftools view -o "${out}" "$@"; then
          return 0
        fi
        echo "$(date -u +%Y-%m-%dT%H:%M:%SZ) bcftools view failed (attempt ${attempt}/${max})" >&2
        if [[ "${attempt}" -eq "${max}" ]]; then
          return 1
        fi
        sleep "${delay}"
        delay=$((delay * 2))
        attempt=$((attempt + 1))
      done
    }

    subset_region() {
      local src="$1"
      local idx="$2"
      local dst="$3"
      local query="$4"
      echo "Subset ${query}: ${src}" >&2
      bcftools_view_retry "${dst}" -r "${query}" -Oz --threads ~{bcftools_threads} \
        "${src}##idx##${idx}"
      bcftools index -t "${dst}"
      echo -n "$(basename "${dst}") records: "
      bcftools index -n "${dst}" 2>/dev/null || bcftools view -H "${dst}" | wc -l
    }

    subset_region "~{flare_vcf}" "~{flare_vcf_index}" flare.region.vcf.gz "${FLARE_REGION}"
    subset_region "~{target_vcf}" "~{target_vcf_index}" target.region.vcf.gz "${REGION}"

    python3 "~{propagate_script}" \
      --flare flare.region.vcf.gz \
      --target target.region.vcf.gz \
      --output annotated.vcf \
      --stats-json "~{prefix}.propagate_stats.json" \
      --stats-tsv "~{prefix}.propagate_stats.tsv" \
      --missing-tsv "~{prefix}.missing.tsv" \
      --tags "~{format_tags}" \
      --min-sample-overlap ~{min_sample_overlap}

    bcftools view -Oz --threads ~{bcftools_threads} -o "~{prefix}.ancestry.vcf.gz" annotated.vcf
    rm -f annotated.vcf flare.region.vcf.gz flare.region.vcf.gz.tbi target.region.vcf.gz target.region.vcf.gz.tbi
    bcftools index -t "~{prefix}.ancestry.vcf.gz"
    test -s "~{prefix}.ancestry.vcf.gz"
    ls -lh "~{prefix}.ancestry.vcf.gz"* "~{prefix}.propagate_stats."* "~{prefix}.missing.tsv" || true
    df -h . || true

    python3 - <<'PY'
import json, sys
from pathlib import Path
stats = json.loads(Path("~{prefix}.propagate_stats.json").read_text())
rate = stats.get("coverage_rate")
thr = float("~{min_coverage_rate}")
n = stats.get("annotated_sites")
n_target = stats.get("target_sites")
exact = stats.get("exact_pos_overlap")
missing = stats.get("sites_missing_ancestry")
print(
    f"coverage_rate={rate} exact_pos={exact} missing={missing} "
    f"({n}/{n_target}) threshold={thr}",
    file=sys.stderr,
)
if rate is not None and thr > 0 and rate < thr:
    msg = f"FLARE ancestry coverage {rate:.4f} < {thr:.4f} ({n}/{n_target} target sites)"
    if "~{fail_on_low_coverage}" == "true":
        raise SystemExit(msg)
    print(f"warning: {msg}; VCF written anyway (set fail_on_low_coverage=true to fail)", file=sys.stderr)
PY
  >>>

  output {
    File annotated_vcf = "~{prefix}.ancestry.vcf.gz"
    File annotated_vcf_index = "~{prefix}.ancestry.vcf.gz.tbi"
    File stats_json = "~{prefix}.propagate_stats.json"
    File stats_tsv = "~{prefix}.propagate_stats.tsv"
    File missing_tsv = "~{prefix}.missing.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
    maxRetries: max_retries
  }
}
