version 1.0

# Per-chromosome signed length × sequence-context spectrum for Figure 2,
# on the DeepVariant + GLnexus joint callset (same shards as Table 2).
#
# Launch from the Terra GL_INTERVAL_set data table (one workflow per row):
#   chrom = this.GL_INTERVAL_set_id
#   vcf   = this.VCF
#
# Counts SNVs (length 0) and indels with |L| < small_max_bp (FILTER PASS or ".",
# as in BcftoolsGlnexusStats). Resolved SVs >= 20 bp, ultralong, and BND come
# from the v3 SV companions, scanned once in MergeVariantLengthSpectrum.
#
# Emits 1 bp signed-length counts per region_class (US / RM / SD / SR, rule of
# annotate_repeat_context.py). Display binning is chosen when plotting.

workflow VariantLengthSpectrum {
  input {
    String chrom
    File vcf

    File rmsk_bed
    File simple_repeat_bed
    File segdup_bed

    File summarize_py
    File sv_site_utils_py
    File annotate_repeat_context_py

    Int small_max_bp = 20
    String apply_filters = "PASS,."

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-sv-annotation:0.1.6"
    Int cpu = 2
    Int memory_gb = 8
    Int preemptible = 2
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 1.2
  }

  call SummarizeChrom {
    input:
      chrom = chrom,
      vcf = vcf,
      rmsk_bed = rmsk_bed,
      simple_repeat_bed = simple_repeat_bed,
      segdup_bed = segdup_bed,
      summarize_py = summarize_py,
      sv_site_utils_py = sv_site_utils_py,
      annotate_repeat_context_py = annotate_repeat_context_py,
      small_max_bp = small_max_bp,
      apply_filters = apply_filters,
      docker = docker,
      cpu = cpu,
      memory_gb = memory_gb,
      preemptible = preemptible,
      disk_gb_floor = disk_gb_floor,
      disk_gb_multiplier = disk_gb_multiplier
  }

  output {
    String chrom_id = SummarizeChrom.chrom_id
    File bins = SummarizeChrom.bins
    File summary = SummarizeChrom.summary
  }

  meta {
    description: "1 bp signed-length × region_class counts from one GLnexus chromosome shard for Figure 2."
    allowNestedInputs: true
  }
}

task SummarizeChrom {
  input {
    String chrom
    File vcf
    File rmsk_bed
    File simple_repeat_bed
    File segdup_bed
    File summarize_py
    File sv_site_utils_py
    File annotate_repeat_context_py
    Int small_max_bp
    String apply_filters
    String docker
    Int cpu
    Int memory_gb
    Int preemptible
    Int disk_gb_floor
    Float disk_gb_multiplier
  }

  Int disk_gb = ceil(size(vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    python3 -c "import pathlib,sys; pathlib.Path('chrom.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "~{chrom}"
    CHROM="$(cat chrom.txt)"
    # Prefer localized scripts over whatever is baked into the image.
    mkdir -p "$PWD/py"
    cp "~{summarize_py}" "$PWD/py/summarize_variant_length_spectrum.py"
    cp "~{sv_site_utils_py}" "$PWD/py/sv_site_utils.py"
    cp "~{annotate_repeat_context_py}" "$PWD/py/annotate_repeat_context.py"
    export PYTHONPATH="$PWD/py:${PYTHONPATH:-}"

    echo "[$(date -Is)] length spectrum ${CHROM}" >&2
    ls -lh "~{vcf}" || true
    df -h . || true

    python3 "$PWD/py/summarize_variant_length_spectrum.py" \
      --glnexus-vcf "~{vcf}" \
      --chrom "${CHROM}" \
      --rmsk-bed "~{rmsk_bed}" \
      --simple-repeat-bed "~{simple_repeat_bed}" \
      --segdup-bed "~{segdup_bed}" \
      --small-max-bp ~{small_max_bp} \
      --apply-filters "~{apply_filters}" \
      --out-bins "${CHROM}.length_spectrum.bins.tsv" \
      --out-summary "${CHROM}.length_spectrum.summary.json"

    ln -s "${CHROM}.length_spectrum.bins.tsv" bins.tsv
    ln -s "${CHROM}.length_spectrum.summary.json" summary.json
  >>>

  output {
    String chrom_id = read_string("chrom.txt")
    File bins = "bins.tsv"
    File summary = "summary.json"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GiB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}
