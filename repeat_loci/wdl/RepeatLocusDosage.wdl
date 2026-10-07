version 1.0

# Per-locus repeat dosage on the integrated SNV/indel/SV callset (Figure 2C).
#
# Sums the signed length change of every record that falls in a catalog
# repeat locus (TRExplorer / Adotto BED from trgt_plvi.py catalog) per
# haplotype, and writes one allele table per chromosome (locus x distinct
# dosage, haplotype counts, no sample ids) plus diagnostics for split and
# double-counted alleles. See scripts/aggregate_repeat_loci.py.
#
# Inputs, either or both:
#   chrom_vcfs / chrom_vcf_indexes   one VCF per entry of `chroms`, localized
#                                    (Phase 2: DeepVariant + GLnexus shards)
#   genome_vcfs / genome_vcf_indexes genome-wide VCF/BCFs, cut per chromosome
#                                    to length-changing records before use
#                                    (Phase 2: v3_main.bcf; Phase 1: SV BCF +
#                                    DeepVariant joint VCF)
#
# Roles (chrom_vcf_role, genome_vcf_roles): "small" keeps alleles with
# |change| < split_bp, "sv" keeps the rest, "all" keeps every allele, so an
# event called by both the small-variant and SV callsets counts once.

workflow RepeatLocusDosage {
  input {
    String label
    Array[String] chroms
    Array[File] chrom_vcfs = []
    Array[File] chrom_vcf_indexes = []
    Array[File] genome_vcfs = []
    Array[File] genome_vcf_indexes = []
    String chrom_vcf_role = "all"
    Array[String] genome_vcf_roles = []
    Int split_bp = 50
    String apply_filters = "PASS,."

    File catalog_bed
    File? samples
    String exclude_prefixes = "HG,NA"
    Int flank = 10
    Int max_bp = 100000
    Boolean with_ancestry = false
    Boolean ignore_phase = false
    Boolean write_locus_vcf = false

    File aggregate_repeat_loci_py

    String gcs_project = ""
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int cpu = 2
    Int memory_gb = 16
    Int preemptible = 1
    Int disk_gb_floor = 50
  }

  Boolean has_chrom_vcfs = length(chrom_vcfs) > 0
  Boolean has_genome_vcfs = length(genome_vcfs) > 0

  scatter (i in range(length(chroms))) {
    String chrom = chroms[i]

    if (has_chrom_vcfs) {
      File chrom_vcf = chrom_vcfs[i]
      File chrom_vcf_index = chrom_vcf_indexes[i]
      String chrom_role = chrom_vcf_role
    }

    if (has_genome_vcfs) {
      scatter (j in range(length(genome_vcfs))) {
        String genome_role = if length(genome_vcf_roles) > 0 then genome_vcf_roles[j] else "all"
        call SliceChrom {
          input:
            vcf = genome_vcfs[j],
            vcf_index = genome_vcf_indexes[j],
            chrom = chrom,
            tag = "g~{j}",
            gcs_project = gcs_project,
            docker = docker,
            preemptible = preemptible
        }
      }
    }

    call AggregateChrom {
      input:
        chrom = chrom,
        label = label,
        vcfs = flatten([select_all([chrom_vcf]), select_first([SliceChrom.bcf, []])]),
        vcf_indexes = flatten([select_all([chrom_vcf_index]), select_first([SliceChrom.bcf_index, []])]),
        vcf_roles = flatten([select_all([chrom_role]), select_first([genome_role, []])]),
        split_bp = split_bp,
        apply_filters = apply_filters,
        catalog_bed = catalog_bed,
        samples = samples,
        exclude_prefixes = exclude_prefixes,
        flank = flank,
        max_bp = max_bp,
        with_ancestry = with_ancestry,
        ignore_phase = ignore_phase,
        write_locus_vcf = write_locus_vcf,
        aggregate_repeat_loci_py = aggregate_repeat_loci_py,
        docker = docker,
        cpu = cpu,
        memory_gb = memory_gb,
        preemptible = preemptible,
        disk_gb_floor = disk_gb_floor
    }
  }

  call MergeChroms {
    input:
      label = label,
      allele_tables = AggregateChrom.alleles,
      summary_files = AggregateChrom.summary,
      docker = docker
  }

  output {
    File alleles = MergeChroms.alleles
    File summary = MergeChroms.summary
    Array[File] chrom_alleles = AggregateChrom.alleles
    Array[File?] locus_vcfs = AggregateChrom.locus_vcf
    Array[File?] locus_vcf_indexes = AggregateChrom.locus_vcf_index
  }

  meta {
    description: "Per-locus repeat dosage alleles from the integrated callset, one task per chromosome."
    allowNestedInputs: true
  }
}

# A Cromwell GCS_OAUTH_TOKEN expires after ~1 h; mint a fresh one at start
# and read the chromosome in one bcftools pass.
task SliceChrom {
  input {
    File vcf
    File vcf_index
    String chrom
    String tag
    String gcs_project
    String docker
    Int preemptible
  }

  parameter_meta {
    vcf: { description: "Genome-wide VCF/BCF; region-streamed, not localized.", localization_optional: true }
    vcf_index: { description: "Tabix/CSI sibling of vcf.", localization_optional: true }
  }

  command <<<
    set -euo pipefail
    export HTS_RETRY_MAX="${HTS_RETRY_MAX:-8}"
    export HTS_RETRY_DELAY="${HTS_RETRY_DELAY:-500}"
    python3 - "~{vcf}" <<'PY'
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
project = "~{gcs_project}".strip() or os.environ.get("GOOGLE_PROJECT", "")
if project:
    parts.append(f"export GCS_REQUESTER_PAYS_PROJECT={project}")
if sys.argv[1].startswith("gs://"):
    token = fresh_token()
    if not token:
        raise SystemExit("error: bcftools needs a GCS OAuth token to read gs:// VCFs")
    parts.append(f"export GCS_OAUTH_TOKEN='{token}'")
Path("gcs.env").write_text(("\n".join(parts) + "\n") if parts else "")
PY
    # shellcheck disable=SC1091
    source gcs.env

    bcftools view -r "~{chrom}" -i 'strlen(REF)!=strlen(ALT) || ALT~"<"' \
      -Ob -o "~{chrom}.~{tag}.bcf" --write-index "~{vcf}##idx##~{vcf_index}"
    ls -lh "~{chrom}.~{tag}.bcf"*
  >>>

  output {
    File bcf = "~{chrom}.~{tag}.bcf"
    File bcf_index = "~{chrom}.~{tag}.bcf.csi"
  }

  runtime {
    docker: docker
    cpu: 2
    memory: "4 GiB"
    disks: "local-disk 50 HDD"
    preemptible: preemptible
    maxRetries: 2
  }
}

task AggregateChrom {
  input {
    String chrom
    String label
    Array[File] vcfs
    Array[File] vcf_indexes
    Array[String] vcf_roles
    Int split_bp
    String apply_filters
    File catalog_bed
    File? samples
    String exclude_prefixes
    Int flank
    Int max_bp
    Boolean with_ancestry
    Boolean ignore_phase
    Boolean write_locus_vcf
    File aggregate_repeat_loci_py
    String docker
    Int cpu
    Int memory_gb
    Int preemptible
    Int disk_gb_floor
  }

  Int disk_gb = ceil(size(vcfs, "GB") * 1.2) + disk_gb_floor
  String prefix = "~{label}.~{chrom}"

  command <<<
    set -euo pipefail
    mkdir -p py in
    cp "~{aggregate_repeat_loci_py}" py/aggregate_repeat_loci.py

    # Pair each VCF with its index under one directory so bcftools finds it.
    paste "~{write_lines(vcfs)}" "~{write_lines(vcf_indexes)}" "~{write_lines(vcf_roles)}" > pairs.tsv
    VCF_ARGS=()
    n=0
    while IFS=$'\t' read -r v x role; do
      n=$((n + 1))
      name="in/${n}.$(basename "${v}")"
      ln -s "${v}" "${name}"
      case "${x}" in *.csi) ln -s "${x}" "${name}.csi" ;; *) ln -s "${x}" "${name}.tbi" ;; esac
      case "${role}" in
        all) VCF_ARGS+=(--vcf "${name}") ;;
        small) VCF_ARGS+=(--small-vcf "${name}") ;;
        sv) VCF_ARGS+=(--sv-vcf "${name}") ;;
        *) echo "unknown VCF role ${role}" >&2; exit 1 ;;
      esac
    done < pairs.tsv
    test "${n}" -gt 0

    echo "[$(date -Is)] repeat loci ~{prefix} (${n} VCFs)" >&2
    python3 py/aggregate_repeat_loci.py \
      "${VCF_ARGS[@]}" \
      --catalog-bed "~{catalog_bed}" \
      --chrom "~{chrom}" \
      ~{"--samples " + samples} \
      --exclude-prefixes "~{exclude_prefixes}" \
      --flank ~{flank} \
      --max-bp ~{max_bp} \
      --split-bp ~{split_bp} \
      --apply-filters "~{apply_filters}" \
      ~{true="--with-ancestry" false="" with_ancestry} \
      ~{true="--ignore-phase" false="" ignore_phase} \
      ~{true="--out-vcf ~{prefix}.loci.vcf" false="" write_locus_vcf} \
      --out-alleles "~{prefix}.alleles.tsv.gz" \
      --out-summary "~{prefix}.summary.json"

    if [[ -f "~{prefix}.loci.vcf" ]]; then
      bgzip -f "~{prefix}.loci.vcf"
      tabix -p vcf "~{prefix}.loci.vcf.gz"
    fi
  >>>

  output {
    File alleles = "~{prefix}.alleles.tsv.gz"
    File summary = "~{prefix}.summary.json"
    File? locus_vcf = "~{prefix}.loci.vcf.gz"
    File? locus_vcf_index = "~{prefix}.loci.vcf.gz.tbi"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task MergeChroms {
  input {
    String label
    Array[File] allele_tables
    Array[File] summary_files
    String docker
  }

  command <<<
    set -euo pipefail
    python3 - "~{label}" "~{write_lines(allele_tables)}" "~{write_lines(summary_files)}" <<'PY'
import gzip, json, sys
from pathlib import Path

label, allele_list, summary_list = sys.argv[1:4]
with gzip.open(f"{label}.alleles.tsv.gz", "wt") as out:
    header = None
    for path in Path(allele_list).read_text().split():
        with gzip.open(path, "rt") as fh:
            h = fh.readline()
            if header is None:
                header = h
                out.write(h)
            elif h != header:
                raise SystemExit(f"header mismatch in {path}")
            for line in fh:
                out.write(line)

counts, per_chrom = {}, {}
for path in Path(summary_list).read_text().split():
    s = json.loads(Path(path).read_text())
    per_chrom[Path(path).name] = s
    for k, v in s["counts"].items():
        counts[k] = counts.get(k, 0) + v
c = counts
merged = {
    "label": label,
    "counts": c,
    "hap_missing_frac": c["hap_missing"] / max(c["hap_called"] + c["hap_missing"], 1),
    "hap_possible_duplicate_frac": c["hap_possible_duplicate"] / max(c["hap_called"], 1),
    "per_chrom": per_chrom,
}
Path(f"{label}.summary.json").write_text(json.dumps(merged, indent=2) + "\n")
PY
  >>>

  output {
    File alleles = "~{label}.alleles.tsv.gz"
    File summary = "~{label}.summary.json"
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "4 GiB"
    disks: "local-disk 50 HDD"
    preemptible: 1
  }
}
