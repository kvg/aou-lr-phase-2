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
    String region

    File bam
    File bam_index
    File phased_vcf
    File phased_vcf_index
    File reference_fasta
    File reference_fasta_index
    File catalog_bed

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
    Int disk_gb = 100
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
      catalog_bed = catalog_bed,
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
    String region
    File bam
    File bam_index
    File phased_vcf
    File phased_vcf_index
    File reference_fasta
    File reference_fasta_index
    File catalog_bed
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
    Int disk_gb
    Int preemptible
  }

  parameter_meta {
    bam: { description: "Haplotagged BAM; only the window is read, through the index.", localization_optional: true }
    bam_index: { description: "BAI/CSI sibling of bam.", localization_optional: true }
    phased_vcf: { description: "Ligated phased callset for the chromosome; read by region.", localization_optional: true }
    phased_vcf_index: { description: "Tabix or CSI sibling of phased_vcf.", localization_optional: true }
  }

  command <<<
    set -euo pipefail
    export HTS_RETRY_MAX="${HTS_RETRY_MAX:-8}"
    export HTS_RETRY_DELAY="${HTS_RETRY_DELAY:-500}"

    python3 - "~{bam}" "~{phased_vcf}" "~{production_trgt_vcf}" "~{gcs_project}" <<'PY'
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
    # shellcheck disable=SC1091
    source gcs.env

    CHROM=$(echo "~{region}" | cut -d: -f1)
    START=$(echo "~{region}" | cut -d: -f2 | cut -d- -f1)
    END=$(echo "~{region}" | cut -d: -f2 | cut -d- -f2)
    PAD_START=$(( START > ~{variant_flank} ? START - ~{variant_flank} : 1 ))
    PAD_END=$(( END + ~{variant_flank} ))
    PADDED="${CHROM}:${PAD_START}-${PAD_END}"
    echo "window ~{region}; variant window ${PADDED}" >&2

    # 1. One sample out of the phased callset. This is the haplotype frame:
    #    slot 0 here is the haplotype AN1 describes. No --force-samples: a
    #    sample that is absent must fail here, not silently phase nothing.
    bcftools view \
      -s "~{sample_id}" \
      -r "${PADDED}" \
      -Oz -o sample.phased.vcf.gz \
      "~{phased_vcf}##idx##~{phased_vcf_index}"
    bcftools index -t sample.phased.vcf.gz
    bcftools index -n sample.phased.vcf.gz > n_phased_sites.txt
    echo "phased sites in window: $(cat n_phased_sites.txt)" >&2

    # 2. Reads for the window, with HiPhase's HP/PS removed. Those tags are in
    #    HiPhase's block frame; whatshap writes the frame we want in step 3.
    samtools view -b -h -x HP,PS \
      "~{bam}##idx##~{bam_index}" "~{region}" > reads.bam
    samtools index reads.bam
    samtools view -c reads.bam > n_reads.txt

    # 3. Haplotag against the phased callset.
    whatshap haplotag \
      --reference "~{reference_fasta}" \
      --sample "~{sample_id}" \
      ~{true="--ignore-read-groups" false="" ignore_read_groups} \
      --output-haplotag-list haplotag.tsv \
      --output tagged.bam \
      sample.phased.vcf.gz \
      reads.bam
    samtools index tagged.bam

    # 4. Catalog rows whose locus starts inside the window. TRGT would accept
    #    the whole catalog, but only the window has reads.
    if [[ "~{catalog_bed}" == *.gz ]]; then DECOMP=zcat; else DECOMP=cat; fi
    ${DECOMP} "~{catalog_bed}" \
      | awk -v c="${CHROM}" -v s="${START}" -v e="${END}" \
          'BEGIN{FS=OFS="\t"} $1==c && $2>=s-1 && $2<=e' > catalog.window.bed
    wc -l < catalog.window.bed > n_catalog.txt
    echo "catalog loci in window: $(cat n_catalog.txt)" >&2
    if [[ ! -s catalog.window.bed ]]; then
      echo "error: no catalog locus starts inside ~{region}" >&2
      exit 1
    fi

    # 5. Re-genotype from the tagged reads. Detect whether this TRGT takes the
    #    `genotype` subcommand (1.x) or the bare form (0.x).
    if trgt --help 2>&1 | grep -qw genotype; then SUB=(genotype); else SUB=(); fi
    trgt "${SUB[@]}" \
      --genome "~{reference_fasta}" \
      --repeats catalog.window.bed \
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
      echo -e "catalog_loci_in_window\t$(cat n_catalog.txt)"
      echo -e "haplotag_rows\t$(( $(wc -l < haplotag.tsv) - 1 ))"
      echo -e "spanning_reads\t$(samtools view -c "~{sample_id}.spanning.bam")"
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

    # The production call for the same window, for the notebook's comparison of
    # re-genotyped against as-delivered calls.
    if [[ -n "~{production_trgt_vcf}" ]]; then
      bcftools view -r "~{region}" -Oz -o "~{sample_id}.production_trgt.vcf.gz" \
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
