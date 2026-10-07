version 1.0

# Checks per-locus repeat dosage from the integrated callset against TRGT on
# one chromosome for a participant subset: exact / within-one-unit agreement
# and integrated-sum overcounts (split or double-counted records), stratified
# by how many integrated records hit the locus. See
# scripts/compare_trgt_locus_dosage.py.
#
# trgt_table: participant<TAB>TRGT VCF URI, the subset to compare.
# vcfs / vcf_roles: the same callsets and roles as RepeatLocusDosage
# (Phase 2: GLnexus shard for `chrom` as "small", v3_main.bcf as "sv").
# region: indexed window. The callsets and TRGT VCFs stay in the bucket;
# bcftools reads this window through the index.

workflow TrgtLocusConcordance {
  input {
    String chrom = "chr20"
    String region = "chr20:10000000-20000000"
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[String] vcf_roles
    File catalog_bed
    File trgt_table
    Int flank = 10
    Int split_bp = 50
    String apply_filters = "PASS,."
    Boolean ignore_phase = true
    String gcs_project = ""

    File aggregate_repeat_loci_py
    File compare_trgt_locus_dosage_py

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int memory_gb = 16
    Int preemptible = 1
  }

  Array[Array[String]] columns = transpose(read_tsv(trgt_table))

  call Concordance {
    input:
      chrom = chrom,
      region = region,
      gcs_project = gcs_project,
      vcfs = vcfs,
      vcf_indexes = vcf_indexes,
      vcf_roles = vcf_roles,
      split_bp = split_bp,
      apply_filters = apply_filters,
      ignore_phase = ignore_phase,
      catalog_bed = catalog_bed,
      samples = columns[0],
      trgt_vcfs = columns[1],
      flank = flank,
      aggregate_repeat_loci_py = aggregate_repeat_loci_py,
      compare_trgt_locus_dosage_py = compare_trgt_locus_dosage_py,
      docker = docker,
      memory_gb = memory_gb,
      preemptible = preemptible
  }

  output {
    File loci = Concordance.loci
    File summary = Concordance.summary
    File aggregate_summary = Concordance.aggregate_summary
  }

  meta {
    description: "Integrated-callset locus dosage vs TRGT on one chromosome (double-count check for Figure 2C)."
    allowNestedInputs: true
  }
}

task Concordance {
  input {
    String chrom
    String region
    String gcs_project
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[String] vcf_roles
    Int split_bp
    String apply_filters
    Boolean ignore_phase
    File catalog_bed
    Array[String] samples
    Array[String] trgt_vcfs
    Int flank
    File aggregate_repeat_loci_py
    File compare_trgt_locus_dosage_py
    String docker
    Int memory_gb
    Int preemptible
  }

  parameter_meta {
    vcfs: { description: "Indexed callsets; the region is streamed, the files are not localized.", localization_optional: true }
    vcf_indexes: { description: "Tabix or CSI sibling of each VCF.", localization_optional: true }
  }

  command <<<
    set -euo pipefail
    export HTS_RETRY_MAX="${HTS_RETRY_MAX:-8}"
    export HTS_RETRY_DELAY="${HTS_RETRY_DELAY:-500}"
    python3 - "~{sep=' ' vcfs}" "~{gcs_project}" <<'PY'
import json, os, subprocess, sys, urllib.request
from pathlib import Path

def fresh_token():
    for cmd in (["gcloud", "auth", "print-access-token"],
                ["gcloud", "auth", "application-default", "print-access-token"]):
        try:
            tok = subprocess.check_output(cmd, text=True, stderr=subprocess.DEVNULL).strip()
            if tok:
                return tok
        except Exception:
            pass
    for url in ("http://metadata.google.internal/computeMetadata/v1/instance/service-accounts/default/token",
                "http://169.254.169.254/computeMetadata/v1/instance/service-accounts/default/token"):
        req = urllib.request.Request(url, headers={"Metadata-Flavor": "Google"})
        try:
            with urllib.request.urlopen(req, timeout=5) as resp:
                tok = json.load(resp).get("access_token") or ""
                if tok:
                    return tok
        except Exception:
            continue
    return ""

parts = []
project = sys.argv[-1].strip() or os.environ.get("GOOGLE_PROJECT", "")
if project:
    parts.append(f"export GCS_REQUESTER_PAYS_PROJECT={project}")
if any(a.startswith("gs://") for a in sys.argv[1:-1]):
    token = fresh_token()
    if not token:
        raise SystemExit("error: bcftools needs a GCS OAuth token to read gs:// VCFs")
    parts.append(f"export GCS_OAUTH_TOKEN='{token}'")
Path("gcs.env").write_text(("\n".join(parts) + "\n") if parts else "")
PY
    # shellcheck disable=SC1091
    source gcs.env

    # Keep the gs:// URI and name its index. bcftools then reads only --region.
    paste "~{write_lines(vcfs)}" "~{write_lines(vcf_indexes)}" "~{write_lines(vcf_roles)}" > pairs.tsv
    VCF_ARGS=()
    while IFS=$'\t' read -r v x role; do
      spec="${v}##idx##${x}"
      case "${role}" in
        all) VCF_ARGS+=(--vcf "${spec}") ;;
        small) VCF_ARGS+=(--small-vcf "${spec}") ;;
        sv) VCF_ARGS+=(--sv-vcf "${spec}") ;;
        *) echo "unknown VCF role ${role}" >&2; exit 1 ;;
      esac
    done < pairs.tsv
    cp "~{write_lines(samples)}" samples.txt
    paste samples.txt "~{write_lines(trgt_vcfs)}" > trgt.tsv

    python3 "~{aggregate_repeat_loci_py}" \
      "${VCF_ARGS[@]}" \
      --catalog-bed "~{catalog_bed}" \
      --region "~{region}" \
      --samples samples.txt \
      --flank ~{flank} \
      --split-bp ~{split_bp} \
      --apply-filters "~{apply_filters}" \
      ~{true="--ignore-phase" false="" ignore_phase} \
      --out-alleles alleles.tsv.gz \
      --out-vcf loci.vcf \
      --out-summary aggregate.summary.json

    python3 "~{compare_trgt_locus_dosage_py}" \
      --locus-vcf loci.vcf \
      --catalog-bed "~{catalog_bed}" \
      --trgt-tsv trgt.tsv \
      --region "~{region}" \
      --out-loci "~{chrom}.trgt_concordance.loci.tsv.gz" \
      --out-summary "~{chrom}.trgt_concordance.summary.json"
  >>>

  output {
    File loci = "~{chrom}.trgt_concordance.loci.tsv.gz"
    File summary = "~{chrom}.trgt_concordance.summary.json"
    File aggregate_summary = "aggregate.summary.json"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: memory_gb + " GiB"
    disks: "local-disk 20 HDD"
    preemptible: preemptible
  }
}
