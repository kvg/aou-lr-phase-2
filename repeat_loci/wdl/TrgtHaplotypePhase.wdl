version 1.0

# Puts TRGT repeat alleles onto the haplotypes of the joint phased callset, for
# one sample over one window, so a per-ancestry repeat dosage can be built on
# TRGT lengths instead of summed SNV/indel/SV records.
#
# Why this exists: the production TRGT VCFs are raw TRGT output. Their
# genotypes are '/'-separated with no PS, so the order of AL / ALLR / SD / MC /
# MS / AP / AM is arbitrary and no allele can be attributed to the haplotype
# carrying a FLARE ancestry call. HiPhase was run on the small-variant and SV
# BCFs only, never on the TRGT VCF, so there is no phased sibling to use.
#
# What it does, and what each step assumes:
#   1. bcftools  pulls the one sample out of the ligated phased callset over
#      the window. This defines the haplotype frame: haplotype 1 here is the
#      haplotype AN1 describes.
#   2. samtools  reads the window out of the haplotagged BAM and STRIPS HP/PS.
#      The production tags come from HiPhase and sit in HiPhase's read-backed
#      block frame, which is not the SHAPEIT4 frame; keeping them would silently
#      mix the two.
#   3. whatshap haplotag re-tags the reads against the per-sample phased VCF.
#      Statistical phase spans the chromosome, so the resulting HP tags are
#      consistent across the window with no block registration to do.
#   4. trgt genotype re-calls the repeats from the tagged reads. TRGT has used
#      HP during genotyping since 0.7.0 and written read-to-allele assignments
#      into the spanning BAM since 0.5.0.
#   5. phase_trgt_by_haplotag.py joins read->haplotype against read->allele and
#      emits the per-haplotype allele table plus a phased TRGT VCF.
#
# This is a prototype: one sample, one window, nothing downstream of it. It
# does not touch FelixGenome.wdl or any cached result.
#
# FIRST RUN: leave allele_tag empty. The task then runs the join in
# --census-only mode, which writes the spanning BAM's SAM-tag census and stops.
# Read the tag that carries the allele index out of phase_summary, then re-run
# with allele_tag set. The tag is deliberately not defaulted: it is
# TRGT-version-dependent, and a wrong guess would produce an empty join rather
# than an error.

workflow TrgtHaplotypePhase {
  input {
    String sample_id
    String region = ""

    File bam
    File bam_index
    File phased_vcf
    File phased_vcf_index
    File reference_fasta
    File reference_fasta_index
    File trgt_repeats_bed

    String production_trgt_vcf = ""
    String production_trgt_vcf_index = ""

    String allele_tag = ""
    Int variant_flank = 20000
    Int min_reads = 2
    Float min_purity = 0.8
    Boolean ignore_read_groups = false
    String gcs_project = ""

    File phase_trgt_by_haplotag_py

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-trgt-phase@sha256:4a6d36d095aa7021431ce65bdea84780be5de58e23683869c01d2b58f890ac88"
    Int cpu = 4
    Int memory_gb = 32
    # Whole-genome mode writes two full copies of the BAM to local disk -- the
    # HP/PS-stripped reads.bam and whatshap's tagged.bam -- plus the spanning
    # BAM and a 3.1 GB reference, so the default is derived from the input
    # rather than guessed. A windowed run reads only its slice and needs none of
    # that headroom.
    Int disk_gb = if region == "" then ceil(2.8 * size(bam, "GB")) + 80 else 100
    Int preemptible = 1
  }

  call HaplotagAndGenotype {
    input:
      sample_id = sample_id,
      region = region,
      bam = bam,
      bam_index = bam_index,
      phased_vcf = phased_vcf,
      phased_vcf_index = phased_vcf_index,
      reference_fasta = reference_fasta,
      reference_fasta_index = reference_fasta_index,
      trgt_repeats_bed = trgt_repeats_bed,
      production_trgt_vcf = production_trgt_vcf,
      production_trgt_vcf_index = production_trgt_vcf_index,
      allele_tag = allele_tag,
      variant_flank = variant_flank,
      min_reads = min_reads,
      min_purity = min_purity,
      ignore_read_groups = ignore_read_groups,
      gcs_project = gcs_project,
      phase_trgt_by_haplotag_py = phase_trgt_by_haplotag_py,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      bam_size_gb = size(bam, "GB"),
      disk_gb = disk_gb,
      preemptible = preemptible
  }

  output {
    File phase_summary = HaplotagAndGenotype.phase_summary
    File? hap_alleles = HaplotagAndGenotype.hap_alleles
    File? phased_trgt_vcf = HaplotagAndGenotype.phased_trgt_vcf
    File trgt_vcf = HaplotagAndGenotype.trgt_vcf
    File trgt_vcf_index = HaplotagAndGenotype.trgt_vcf_index
    File spanning_bam = HaplotagAndGenotype.spanning_bam
    File haplotag_list = HaplotagAndGenotype.haplotag_list
    File read_accounting = HaplotagAndGenotype.read_accounting
    File? production_trgt_slice = HaplotagAndGenotype.production_trgt_slice
  }

  meta {
    description: "Haplotag reads against the phased callset, re-genotype repeats with TRGT, and put each TRGT allele on a phased haplotype (one sample, one window)."
    allowNestedInputs: true
  }
}

task HaplotagAndGenotype {
  input {
    String sample_id
    String region = ""
    File bam
    File bam_index
    File phased_vcf
    File phased_vcf_index
    File reference_fasta
    File reference_fasta_index
    File trgt_repeats_bed
    String production_trgt_vcf
    String production_trgt_vcf_index
    String allele_tag
    Int variant_flank
    Int min_reads
    Float min_purity
    Boolean ignore_read_groups
    String gcs_project
    File phase_trgt_by_haplotag_py
    String docker
    Int cpu
    Int memory_gb
    Float bam_size_gb = 0.0
    Int disk_gb
    Int preemptible
  }

  parameter_meta {
    bam: { description: "Haplotagged BAM; streamed through the index -- only the window is read, or the whole file when region is empty.", localization_optional: true }
    bam_index: { description: "BAI/CSI sibling of bam.", localization_optional: true }
    phased_vcf: { description: "Phased callset carrying this sample. Streamed by region. For the cohort run give the PER-SAMPLE split, not the multi-sample callset: a region read of a multi-sample file pays for all 12,261 samples' genotypes to obtain one.", localization_optional: true }
    phased_vcf_index: { description: "Tabix or CSI sibling of phased_vcf.", localization_optional: true }
  }

  command <<<
    set -euo pipefail
    export HTS_RETRY_MAX="${HTS_RETRY_MAX:-8}"
    export HTS_RETRY_DELAY="${HTS_RETRY_DELAY:-500}"

    # GCS access tokens last about an hour. Whole-genome mode runs far longer
    # than that, and the production-slice step at the very end reads gs://
    # through htslib -- which failed with "Operation not permitted" once the
    # startup token went stale. So the minting runs from a script that can be
    # re-invoked, and every gs:// read refreshes immediately beforehand rather
    # than trusting a token minted at task start.
    cat > mint_gcs_env.py <<'PY'
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
        raise SystemExit("error: htslib needs a GCS OAuth token to read gs:// inputs")
    parts.append(f"export GCS_OAUTH_TOKEN='{token}'")
Path("gcs.env").write_text(("\n".join(parts) + "\n") if parts else "")
PY

    refresh_gcs() {
      python3 mint_gcs_env.py "~{bam}" "~{phased_vcf}" "~{production_trgt_vcf}" "~{gcs_project}"
      # shellcheck disable=SC1091
      source gcs.env
    }
    refresh_gcs

    # An empty region means the whole genome: no windowing anywhere below. The
    # cohort run uses that mode, because TR loci sit one per ~698 bp and so a
    # "only the loci we need" read IS a whole-genome read -- there is nothing to
    # be saved by cutting it into pieces, and one task per sample streams the
    # BAM once rather than once per shard.
    REGION="~{region}"
    if [[ -z "${REGION}" ]]; then
      WHOLE_GENOME=1
      echo "no region given: whole-genome mode" >&2
      # Fail now, not an hour in. Whole-genome mode writes reads.bam and
      # tagged.bam, each about the size of the source BAM, so ENOSPC midway
      # through costs the entire task rather than a few seconds here.
      BAM_GB=$(printf "%.0f" "~{bam_size_gb}")
      FREE_GB=$(df -BG --output=avail . 2>/dev/null | tail -1 | tr -dc '0-9' || echo 0)
      echo "bam ${BAM_GB} GiB; free disk ${FREE_GB} GiB" >&2
      if [[ "${BAM_GB}" -le 0 ]]; then
        echo "warning: BAM size came back 0, so the disk precheck is skipped." >&2
        echo "         size() on a localization_optional input can do this; if the" >&2
        echo "         task later dies on ENOSPC, set disk_gb explicitly." >&2
      fi
      if [[ "${BAM_GB}" -gt 0 && "${FREE_GB}" -gt 0 ]]; then
        NEED_GB=$(( BAM_GB * 2 + 20 ))
        if [[ "${FREE_GB}" -lt "${NEED_GB}" ]]; then
          echo "error: whole-genome mode needs about ${NEED_GB} GiB of local disk" >&2
          echo "       (two BAM-sized copies plus the reference) and only" >&2
          echo "       ${FREE_GB} GiB is free. Raise disk_gb." >&2
          exit 1
        fi
      fi
    else
      WHOLE_GENOME=0
      case "${REGION}" in
        *:*-*) ;;
        *)
          echo "error: region must be 'chrom:start-end' or empty for whole genome." >&2
          echo "       Got '${REGION}'. A bare contig is rejected rather than" >&2
          echo "       silently parsed as start=end='', which would window to nothing." >&2
          exit 1 ;;
      esac
      CHROM=$(echo "${REGION}" | cut -d: -f1)
      START=$(echo "${REGION}" | cut -d: -f2 | cut -d- -f1)
      END=$(echo "${REGION}" | cut -d: -f2 | cut -d- -f2)
      PAD_START=$(( START > ~{variant_flank} ? START - ~{variant_flank} : 1 ))
      PAD_END=$(( END + ~{variant_flank} ))
      PADDED="${CHROM}:${PAD_START}-${PAD_END}"
      echo "window ${REGION}; variant window ${PADDED}" >&2
    fi

    # 1. One sample out of the phased callset. This is the haplotype frame:
    #    slot 0 here is the haplotype AN1 describes. No --force-samples: a
    #    sample that is absent must fail here, not silently phase nothing.
    if [[ "${WHOLE_GENOME}" == "1" ]]; then REGION_ARGS=(); else REGION_ARGS=(-r "${PADDED}"); fi
    refresh_gcs
    bcftools view \
      -s "~{sample_id}" \
      "${REGION_ARGS[@]}" \
      -Oz -o sample.phased.vcf.gz \
      "~{phased_vcf}##idx##~{phased_vcf_index}"
    bcftools index -t sample.phased.vcf.gz
    bcftools index -n sample.phased.vcf.gz > n_phased_sites.txt
    echo "phased sites in window: $(cat n_phased_sites.txt)" >&2

    # 2. Reads for the window, with HiPhase's HP/PS removed. Those tags are in
    #    HiPhase's block frame; whatshap writes the frame we want in step 3.
    if [[ "${WHOLE_GENOME}" == "1" ]]; then READ_REGION=(); else READ_REGION=("${REGION}"); fi
    refresh_gcs
    samtools view -b -h -x HP,PS -@ ~{cpu} \
      "~{bam}##idx##~{bam_index}" "${READ_REGION[@]}" > reads.bam
    samtools index reads.bam
    samtools view -c reads.bam > n_reads.txt

    # 3. Haplotag against the phased callset.
    #
    #    whatshap REFUSES outright when the BAM carries reads on a contig the VCF
    #    header does not declare -- it is not a warning and the untagged reads are
    #    not simply skipped. That bites whenever the phased VCF covers fewer
    #    contigs than the BAM: a single-chromosome callset against a whole-genome
    #    BAM, or a per-sample split with no variants on some contig.
    #
    #    --skip-missing-contigs is whatshap's own bypass and gives the semantics
    #    we want (those reads stay untagged). But it would also mask handing over
    #    the wrong VCF entirely, so the mismatch is MEASURED first and written to
    #    the accounting file. A whole-genome run against a chr20-only callset then
    #    shows ~98% of reads on absent contigs instead of quietly phasing 2%.
    bcftools view -h sample.phased.vcf.gz \
      | sed -n 's/^##contig=<ID=\([^,>]*\).*/\1/p' \
      | sort -u > vcf_contigs.txt
    samtools idxstats reads.bam | awk -F'\t' '$1 != "*" && $3 + $4 > 0 {print $1}' \
      | sort -u > bam_contigs.txt
    comm -23 bam_contigs.txt vcf_contigs.txt > contigs_missing_from_vcf.txt
    N_MISSING=$(wc -l < contigs_missing_from_vcf.txt)
    READS_MISSING=$(samtools idxstats reads.bam \
      | awk -F'\t' 'NR==FNR {miss[$1]=1; next} $1 != "*" && ($1 in miss) {n += $3 + $4}
                     END {print n + 0}' contigs_missing_from_vcf.txt -)
    echo "contigs with reads: $(wc -l < bam_contigs.txt); in VCF header: $(wc -l < vcf_contigs.txt);" \
         "absent from VCF: ${N_MISSING} carrying ${READS_MISSING} reads" >&2
    if [[ "${N_MISSING}" -gt 0 ]]; then
      echo "note: those reads cannot be haplotagged and will be left untagged:" >&2
      sed 's/^/        /' contigs_missing_from_vcf.txt | head -5 >&2
    fi

    whatshap haplotag \
      --reference "~{reference_fasta}" \
      --sample "~{sample_id}" \
      --skip-missing-contigs \
      ~{true="--ignore-read-groups" false="" ignore_read_groups} \
      --output-haplotag-list haplotag.tsv \
      --output tagged.bam \
      sample.phased.vcf.gz \
      reads.bam
    samtools index tagged.bam

    # 4. Repeat definitions whose locus starts inside the window. TRGT would
    #    accept the whole file, but only the window has reads.
    #
    #    This must be TRGT's OWN definition format -- 'chrom start end info',
    #    four tab-separated fields, info being 'ID=...;MOTIFS=...;STRUC=...'.
    #    It is NOT refs/repeat_loci/trexplorer_v1.0.1.catalog.bed.gz, which is a
    #    FIVE-column file ('chrom start end TRID motifs') that trgt_plvi.py
    #    derived from a TRGT VCF for aggregate_repeat_loci.py and
    #    trgt_allele_counts.py. That file has dropped STRUC, which cannot be
    #    reconstructed from a motif list at a compound locus, so it cannot be
    #    converted here. Use the definitions the production TRGT run genotyped
    #    against: that is also what makes these calls locus-for-locus
    #    comparable with the delivered ones.
    if [[ "~{trgt_repeats_bed}" == *.gz ]]; then DECOMP=zcat; else DECOMP=cat; fi
    ${DECOMP} "~{trgt_repeats_bed}" | grep -v '^#' | head -1 > first_def.bed || true
    NF_DEF=$(awk -F'\t' '{print NF; exit}' first_def.bed)
    if [[ "${NF_DEF}" != "4" ]] || ! grep -q 'ID=' first_def.bed; then
      echo "error: --repeats needs TRGT's definition format: 4 tab-separated fields," >&2
      echo "       'chrom start end ID=<trid>;MOTIFS=<m>;STRUC=<s>'." >&2
      echo "       Got ${NF_DEF} field(s) on the first row:" >&2
      sed 's/^/         /' first_def.bed >&2
      echo "       A 5-field row means this is the derived catalog" >&2
      echo "       (chrom start end TRID motifs) that trgt_plvi.py writes for" >&2
      echo "       aggregate_repeat_loci.py. It has no STRUC and cannot be used here." >&2
      exit 1
    fi

    if [[ "${WHOLE_GENOME}" == "1" ]]; then
      ${DECOMP} "~{trgt_repeats_bed}" | grep -v '^#' > repeats.window.bed
    else
      ${DECOMP} "~{trgt_repeats_bed}" \
        | awk -v c="${CHROM}" -v s="${START}" -v e="${END}" \
            'BEGIN{FS=OFS="\t"} /^#/ {next} $1==c && $2>=s-1 && $2<=e' > repeats.window.bed
    fi
    wc -l < repeats.window.bed > n_catalog.txt
    echo "repeat definitions used: $(cat n_catalog.txt)" >&2
    if [[ ! -s repeats.window.bed ]]; then
      echo "error: no repeat definition selected for '${REGION:-<whole genome>}'" >&2
      exit 1
    fi

    # 5. Re-genotype from the tagged reads. Detect whether this TRGT takes the
    #    `genotype` subcommand (1.x) or the bare form (0.x).
    if trgt --help 2>&1 | grep -qw genotype; then SUB=(genotype); else SUB=(); fi
    trgt "${SUB[@]}" \
      --genome "~{reference_fasta}" \
      --repeats repeats.window.bed \
      --reads tagged.bam \
      --threads ~{cpu} \
      --output-prefix trgt_out

    bcftools sort -Oz -o "~{sample_id}.trgt.vcf.gz" trgt_out.vcf.gz
    bcftools index -t "~{sample_id}.trgt.vcf.gz"
    samtools sort -@ ~{cpu} -o "~{sample_id}.spanning.bam" trgt_out.spanning.bam
    samtools index "~{sample_id}.spanning.bam"

    # Read accounting, so an empty join can be told from an untagged one.
    {
      echo -e "metric\tvalue"
      echo -e "reads_in_window\t$(cat n_reads.txt)"
      echo -e "phased_sites_in_variant_window\t$(cat n_phased_sites.txt)"
      echo -e "repeat_definitions_in_window\t$(cat n_catalog.txt)"
      echo -e "haplotag_rows\t$(( $(wc -l < haplotag.tsv) - 1 ))"
      echo -e "spanning_reads\t$(samtools view -c "~{sample_id}.spanning.bam")"
      echo -e "contigs_with_reads\t$(wc -l < bam_contigs.txt)"
      echo -e "contigs_in_vcf_header\t$(wc -l < vcf_contigs.txt)"
      echo -e "contigs_absent_from_vcf\t${N_MISSING}"
      echo -e "reads_on_absent_contigs\t${READS_MISSING}"
    } > "~{sample_id}.read_accounting.tsv"

    # 6. Join. With no allele_tag this writes the tag census and stops, which
    #    is how the tag is discovered on the first run.
    if [[ -z "~{allele_tag}" ]]; then
      echo "allele_tag empty: running the join in --census-only mode" >&2
      python3 "~{phase_trgt_by_haplotag_py}" \
        --trgt-vcf "~{sample_id}.trgt.vcf.gz" \
        --spanning-bam "~{sample_id}.spanning.bam" \
        --census-records 0 \
        --census-only \
        --out-summary "~{sample_id}.phase_summary.json"
    else
      python3 "~{phase_trgt_by_haplotag_py}" \
        --trgt-vcf "~{sample_id}.trgt.vcf.gz" \
        --spanning-bam "~{sample_id}.spanning.bam" \
        --haplotag-list haplotag.tsv \
        --allele-tag "~{allele_tag}" \
        --min-reads ~{min_reads} \
        --min-purity ~{min_purity} \
        --out-alleles "~{sample_id}.hap_alleles.tsv.gz" \
        --out-vcf "~{sample_id}.trgt.phased.vcf" \
        --out-summary "~{sample_id}.phase_summary.json"
      bgzip -f "~{sample_id}.trgt.phased.vcf"
      bcftools index -t "~{sample_id}.trgt.phased.vcf.gz" || true
    fi

    gzip -f haplotag.tsv

    # The production call over the same span, for the notebook's comparison of
    # re-genotyped against as-delivered calls. Last step in the task, so by now
    # the startup token is long dead -- refresh before reading, and pass no -r
    # at all in whole-genome mode (an empty -r "" selects nothing).
    if [[ -n "~{production_trgt_vcf}" ]]; then
      refresh_gcs
      bcftools view "${REGION_ARGS[@]}" -Oz -o "~{sample_id}.production_trgt.vcf.gz" \
        "~{production_trgt_vcf}##idx##~{production_trgt_vcf_index}"
      bcftools index -t "~{sample_id}.production_trgt.vcf.gz"
    fi
  >>>

  output {
    File phase_summary = "~{sample_id}.phase_summary.json"
    File? hap_alleles = "~{sample_id}.hap_alleles.tsv.gz"
    File? phased_trgt_vcf = "~{sample_id}.trgt.phased.vcf.gz"
    File trgt_vcf = "~{sample_id}.trgt.vcf.gz"
    File trgt_vcf_index = "~{sample_id}.trgt.vcf.gz.tbi"
    File spanning_bam = "~{sample_id}.spanning.bam"
    File haplotag_list = "haplotag.tsv.gz"
    File read_accounting = "~{sample_id}.read_accounting.tsv"
    File? production_trgt_slice = "~{sample_id}.production_trgt.vcf.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " SSD"
    preemptible: preemptible
  }
}
