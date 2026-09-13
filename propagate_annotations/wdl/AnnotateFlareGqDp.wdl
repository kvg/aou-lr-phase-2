version 1.0

# Copy FORMAT/GQ, FORMAT/DP, and FORMAT/RNC from a DeepVariant/GLnexus VCF onto FLARE sites.
#
# FLARE is the backbone (GT + AN1/AN2 kept). Both inputs are tabix/csi-indexed.
# The task subsets FLARE with bcftools view -r, restricts DeepVariant to those
# positions (one -r range-read, then -T in the pipe), and streams the annotator
# into bcftools view -Oz.
#
# On Terra (Pipelines API), VCF/index File inputs are localization_optional:
# Cromwell passes gs:// URIs and bcftools/htslib range-reads the region. The
# command also works with local paths (miniwdl / if Cromwell does localize).
# See https://cromwell.readthedocs.io/en/latest/optimizations/FileLocalization/
#
# Designed for a 1 Mb smoke test of ancestry-switch vs genotype-quality QC.
# Pass a full contig (e.g. chr22) to annotate a whole chromosome.
#
# Disk: FLARE extract + DeepVariant sites at FLARE positions + bgzip output.
# The annotator writes VCF text to stdout; bcftools bgzips in the pipe so the
# uncompressed VCF never lands on disk (that used to fill a 40 GB VM).
#
#   flare_vcf            = this.model_chr_anc_vcf   (aou_lr_chrom)
#   deepvariant_vcf      = this.VCF                 (GL_INTERVAL_set)
#   annotate_script      = gs://$WORKSPACE_BUCKET/scripts/annotate_flare_gq_dp.py
#   gcs_project          = Terra google/billing project (requester-pays)

workflow AnnotateFlareGqDp {
  input {
    File flare_vcf
    File flare_vcf_index
    File deepvariant_vcf
    File deepvariant_vcf_index
    String region
    String prefix = "flare_gq_dp"
    File annotate_script
    String gcs_project = ""

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Float min_match_rate = 0.95
    Boolean fail_on_low_match = false
    Float min_sample_overlap = 0.95
    String format_tags = "GQ,DP,RNC"
    Int cpu = 4
    Int memory_gb = 8
    Int preemptible = 1
    # Local FLARE extract + FLARE-position DV extract + bgzip output.
    # Multiplier covers full-contig extracts (and miniwdl localizing inputs).
    Int disk_gb_floor = 80
    Float disk_gb_multiplier = 1.0
    Int bcftools_threads = 4
  }

  call AnnotateRegion {
    input:
      flare_vcf = flare_vcf,
      flare_vcf_index = flare_vcf_index,
      deepvariant_vcf = deepvariant_vcf,
      deepvariant_vcf_index = deepvariant_vcf_index,
      region = region,
      prefix = prefix,
      annotate_script = annotate_script,
      gcs_project = gcs_project,
      docker = docker,
      min_match_rate = min_match_rate,
      fail_on_low_match = fail_on_low_match,
      min_sample_overlap = min_sample_overlap,
      format_tags = format_tags,
      cpu = cpu,
      memory_gb = memory_gb,
      preemptible = preemptible,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier,
      bcftools_threads = bcftools_threads
  }

  output {
    File annotated_vcf = AnnotateRegion.annotated_vcf
    File annotated_vcf_index = AnnotateRegion.annotated_vcf_index
    File stats_json = AnnotateRegion.stats_json
    File stats_tsv = AnnotateRegion.stats_tsv
    File unmatched_tsv = AnnotateRegion.unmatched_tsv
  }

  meta {
    description: "Region-subset FLARE VCF annotated with DeepVariant GQ/DP/RNC (linear streaming merge)."
    allowNestedInputs: true
  }
}

task AnnotateRegion {
  input {
    File flare_vcf
    File flare_vcf_index
    File deepvariant_vcf
    File deepvariant_vcf_index
    String region
    String prefix
    File annotate_script
    String gcs_project
    String docker
    Float min_match_rate
    Boolean fail_on_low_match
    Float min_sample_overlap
    String format_tags
    Int cpu
    Int memory_gb
    Int preemptible
    Int disk_gb_floor
    Float disk_gb_multiplier
    Int bcftools_threads
  }

  Int disk_gb = ceil((size(flare_vcf, "GB") + size(deepvariant_vcf, "GB")) * disk_gb_multiplier) + disk_gb_floor

  parameter_meta {
    flare_vcf: {
      description: "FLARE ancestry VCF. Streamed from GCS on Terra PAPI; must also work as a local path.",
      localization_optional: true
    }
    flare_vcf_index: {
      description: "Tabix/CSI sibling of flare_vcf.",
      localization_optional: true
    }
    deepvariant_vcf: {
      description: "DeepVariant/GLnexus joint VCF with FORMAT/GQ, FORMAT/DP, and FORMAT/RNC.",
      localization_optional: true
    }
    deepvariant_vcf_index: {
      description: "Tabix/CSI sibling of deepvariant_vcf.",
      localization_optional: true
    }
  }

  command <<<
    set -euo pipefail
    python3 -c "import pathlib,sys; pathlib.Path('region.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{region}"
    REGION="$(cat region.txt)"
    echo "[$(date -Is)] AnnotateFlareGqDp ${REGION}" >&2
    echo "FLARE URI: ~{flare_vcf}" >&2
    echo "DV URI: ~{deepvariant_vcf}" >&2
    df -h . || true
    bcftools --version | head -n 2 >&2 || true

    ensure_gcs_auth() {
      local uri="$1"
      case "${uri}" in
        gs://*) ;;
        *) return 0 ;;
      esac
      if [[ -z "${GCS_OAUTH_TOKEN:-}" ]]; then
        if command -v gcloud >/dev/null 2>&1; then
          GCS_OAUTH_TOKEN="$(gcloud auth application-default print-access-token 2>/dev/null || gcloud auth print-access-token 2>/dev/null || true)"
        fi
      fi
      if [[ -z "${GCS_OAUTH_TOKEN:-}" ]]; then
        GCS_OAUTH_TOKEN="$(python3 - <<'PY'
import json, urllib.request
urls = (
    "http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
    "http://169.254.169.254/computeMetadata/v1/instance/service-accounts/default/token",
)
for url in urls:
    req = urllib.request.Request(url, headers={"Metadata-Flavor": "Google"})
    try:
        with urllib.request.urlopen(req, timeout=5) as resp:
            print(json.load(resp)["access_token"])
            break
    except Exception:
        continue
PY
)"
      fi
      if [[ -z "${GCS_OAUTH_TOKEN:-}" ]]; then
        echo "error: bcftools needs GCS_OAUTH_TOKEN to read ${uri}" >&2
        exit 1
      fi
      export GCS_OAUTH_TOKEN
    }

    PROJECT="~{gcs_project}"
    if [[ -z "${PROJECT}" && -n "${GOOGLE_PROJECT:-}" ]]; then
      PROJECT="${GOOGLE_PROJECT}"
    fi
    if [[ -n "${PROJECT}" ]]; then
      export GCS_REQUESTER_PAYS_PROJECT="${PROJECT}"
      echo "GCS_REQUESTER_PAYS_PROJECT=${PROJECT}" >&2
    fi

    ensure_gcs_auth "~{flare_vcf}"
    ensure_gcs_auth "~{deepvariant_vcf}"

    subset_region() {
      local src="$1"
      local idx="$2"
      local dst="$3"
      echo "Subset ${REGION}: ${src}" >&2
      bcftools view -r "${REGION}" -Oz --threads ~{bcftools_threads} \
        -o "${dst}" "${src}##idx##${idx}"
      bcftools index -t "${dst}"
      echo -n "$(basename "${dst}") records: "
      bcftools index -n "${dst}" 2>/dev/null || bcftools view -H "${dst}" | wc -l
    }

    subset_region "~{flare_vcf}" "~{flare_vcf_index}" flare.region.vcf.gz
    ls -lh flare.region.vcf.gz* >&2 || true
    df -h . || true

    # GLnexus chrom shards dwarf FLARE. Keep DV records at FLARE positions
    # only: one -r range-read from GCS, then -T in the pipe (no full DV copy).
    bcftools query -f '%CHROM\t%POS\n' flare.region.vcf.gz | uniq > flare.sites.tsv
    echo "FLARE target positions: $(wc -l < flare.sites.tsv)" >&2
    echo "Subset ${REGION} DeepVariant at FLARE positions: ~{deepvariant_vcf}" >&2
    if [[ ! -s flare.sites.tsv ]]; then
      echo "No FLARE sites in ${REGION}; writing header-only DeepVariant extract" >&2
      bcftools view -h "~{deepvariant_vcf}##idx##~{deepvariant_vcf_index}" \
        | bcftools view -Oz --threads ~{bcftools_threads} -o dv.region.vcf.gz
    else
      bcftools view -r "${REGION}" -Ou --threads ~{bcftools_threads} \
        "~{deepvariant_vcf}##idx##~{deepvariant_vcf_index}" \
        | bcftools view -T flare.sites.tsv -Oz --threads ~{bcftools_threads} \
          -o dv.region.vcf.gz
    fi
    bcftools index -t dv.region.vcf.gz
    echo -n "dv.region.vcf.gz records: "
    bcftools index -n dv.region.vcf.gz 2>/dev/null || bcftools view -H dv.region.vcf.gz | wc -l
    ls -lh dv.region.vcf.gz* flare.sites.tsv >&2 || true
    df -h . || true

    python3 "~{annotate_script}" \
      --flare flare.region.vcf.gz \
      --deepvariant dv.region.vcf.gz \
      --output - \
      --stats-json "~{prefix}.annotate_stats.json" \
      --stats-tsv "~{prefix}.annotate_stats.tsv" \
      --unmatched-tsv "~{prefix}.unmatched.tsv" \
      --tags "~{format_tags}" \
      --min-sample-overlap ~{min_sample_overlap} \
      | bcftools view -Oz --threads ~{bcftools_threads} -o "~{prefix}.flare.gq_dp.vcf.gz"
    rm -f flare.region.vcf.gz flare.region.vcf.gz.tbi dv.region.vcf.gz dv.region.vcf.gz.tbi flare.sites.tsv
    bcftools index -t "~{prefix}.flare.gq_dp.vcf.gz"
    test -s "~{prefix}.flare.gq_dp.vcf.gz"
    ls -lh "~{prefix}.flare.gq_dp.vcf.gz"* "~{prefix}.annotate_stats."* "~{prefix}.unmatched.tsv" || true
    df -h . || true

    python3 - <<'PY'
import json, sys
from pathlib import Path
stats = json.loads(Path("~{prefix}.annotate_stats.json").read_text())
rate = stats.get("match_rate")
thr = float("~{min_match_rate}")
n = stats.get("matched_sites")
n_flare = stats.get("flare_sites")
exact = stats.get("matched_exact")
multi = stats.get("matched_multiallelic")
unmatched = stats.get("unmatched_flare_sites")
print(
    f"match_rate={rate} exact={exact} multiallelic={multi} "
    f"unmatched={unmatched} ({n}/{n_flare}) threshold={thr}",
    file=sys.stderr,
)
if rate is not None and thr > 0 and rate < thr:
    msg = f"FLARE–DeepVariant site match rate {rate:.4f} < {thr:.4f} ({n}/{n_flare} sites)"
    if "~{fail_on_low_match}" == "true":
        raise SystemExit(msg)
    print(f"warning: {msg}; VCF written anyway (set fail_on_low_match=true to fail)", file=sys.stderr)
PY
  >>>

  output {
    File annotated_vcf = "~{prefix}.flare.gq_dp.vcf.gz"
    File annotated_vcf_index = "~{prefix}.flare.gq_dp.vcf.gz.tbi"
    File stats_json = "~{prefix}.annotate_stats.json"
    File stats_tsv = "~{prefix}.annotate_stats.tsv"
    File unmatched_tsv = "~{prefix}.unmatched.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }
}
