version 1.0

# TEMPORARY resume workflow for FelixPilot: Step2 -> concat -> summarize, starting from the
# outputs of an earlier FelixPilot run whose MakeGRM, Pack and Null tasks already succeeded.
# It exists so a Step2 fix can be tested without redoing the ~6 h GRM build.
#
# RunFelixStep2, ConcatPhenotypeScores and SummarizeFelixResults below are copied verbatim from
# FelixPilot.wdl, so results match what a full run would produce. Build the inputs JSON with
# felix/scripts/make_resume_inputs.py. Delete this file once a full FelixPilot run is cached.

task RunFelixStep2 {
  input {
    File packed_tar
    String chrom
    String phenotype
    File null_rda
    File variance_ratio
    File samples_used
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    File run_step2_script
    Int num_ancs = 5
    Int min_mac = 50
    Float pvalcutoff_of_haplotype = 0.05
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 16
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 4.0
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(packed_tar, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p packed step2
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    tar xzf "~{packed_tar}" -C packed
    PREFIX="packed/${CHROM}"
    if [[ ! -f "${PREFIX}.meta" ]]; then
      echo "missing ${PREFIX}.meta after unpack" >&2
      ls -lh packed/ || true
      exit 1
    fi

    ln -sf "~{null_rda}" "step2/${PHENO}.rda"
    ln -sf "~{variance_ratio}" "step2/${PHENO}.varianceRatio.txt"

    Rscript "~{run_step2_script}" \
      --felixla-prefix "${PREFIX}" \
      --chrom "${CHROM}" \
      --null-prefix "step2/${PHENO}" \
      --sample-file "~{samples_used}" \
      --sparse-grm "~{sparse_grm_mtx}" \
      --sparse-grm-ids "~{sparse_grm_sample_ids}" \
      --min-mac ~{min_mac} \
      --n-ancestries ~{num_ancs} \
      --pvalcutoff-of-haplotype ~{pvalcutoff_of_haplotype} \
      --n-threads 1 \
      --out-tsv "${PHENO}.felix.tsv" \
      --out-raw "${PHENO}.felix.raw.txt"

    cp "${PHENO}.felix.tsv" results.felix.tsv
    printf '%s\n' "${PHENO}.felix.tsv" > results.output_name.txt
  >>>

  output {
    File results_tsv = "results.felix.tsv"
    File results_raw = "~{phenotype}.felix.raw.txt"
    File results_name = "results.output_name.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task ConcatPhenotypeScores {
  input {
    String phenotype
    Array[File] shard_tsvs
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 1
    Int memory_gb = 4
    Int disk_gb_floor = 20
    Float disk_gb_multiplier = 2.5
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(shard_tsvs, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    OUT="${PHENO}.felix.tsv"
    export PHENO OUT
    LIST=shards.txt
    cat > "${LIST}" <<'EOF'
~{sep="\n" shard_tsvs}
EOF

    python3 - <<'PY'
from pathlib import Path
import os

pheno = os.environ["PHENO"]
out = Path(os.environ["OUT"])
shards = [Path(l.strip()) for l in Path("shards.txt").read_text().splitlines() if l.strip()]
if not shards:
    raise SystemExit("no shard TSVs to concatenate")

header = None
nrows = 0
tmp = Path(str(out) + ".unsorted")
with tmp.open("w", encoding="utf-8") as fo:
    for path in shards:
        with path.open(encoding="utf-8") as fi:
            h = fi.readline()
            if not h:
                raise SystemExit(f"empty shard: {path}")
            if header is None:
                header = h
                fo.write(header)
            elif h != header:
                raise SystemExit(
                    f"header mismatch in {path}:\n  got: {h!r}\n  exp: {header!r}"
                )
            for line in fi:
                if line.strip():
                    fo.write(line)
                    nrows += 1

hdr = header.rstrip("\n")
body = tmp.read_text(encoding="utf-8").splitlines()[1:]

def _sort_key(line: str):
    parts = line.split("\t", 2)
    chrom = parts[0] if parts else ""
    pos = parts[1] if len(parts) > 1 else ""
    return (chrom, int(pos) if pos.isdigit() else pos)

body_sorted = sorted(body, key=_sort_key)
with out.open("w", encoding="utf-8") as fo:
    fo.write(hdr + "\n")
    for line in body_sorted:
        fo.write(line + "\n")
tmp.unlink()
print(f"Wrote {out} ({nrows} variant rows from {len(shards)} shards for {pheno})")
PY
    wc -l "${OUT}"
    cp "${OUT}" merged.felix.tsv
    printf '%s\n' "${OUT}" > merged.output_name.txt
  >>>

  output {
    File merged_tsv = "merged.felix.tsv"
    File merged_name = "merged.output_name.txt"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task SummarizeFelixResults {
  input {
    Array[File] result_tsvs
    Array[String] phenotypes
    File pheno_cov
    File summarize_script
    Float p_threshold = 0.00000005
    Int top_n = 50
    String p_column = "P_cct_admixed_c"
    String named_suffix = ".felix.tsv"
    String report_title = "FELIX chr22 QC summary"
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 2
    Int memory_gb = 16
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 2.0
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(result_tsvs, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p summary
    cat > phenotypes.txt <<'EOF'
~{sep="\n" phenotypes}
EOF
    cat > results.txt <<'EOF'
~{sep="\n" result_tsvs}
EOF

    mapfile -t RESULT_ARR < results.txt
    mapfile -t PHENO_ARR < phenotypes.txt
    if [[ "${#RESULT_ARR[@]}" -ne "${#PHENO_ARR[@]}" ]]; then
      echo "results (${#RESULT_ARR[@]}) != phenotypes (${#PHENO_ARR[@]})" >&2
      exit 1
    fi

    python3 "~{summarize_script}" \
      --results "${RESULT_ARR[@]}" \
      --phenotypes "${PHENO_ARR[@]}" \
      --pheno-cov "~{pheno_cov}" \
      --out-dir summary \
      --top-n ~{top_n} \
      --p-threshold ~{p_threshold} \
      --p-column "~{p_column}" \
      --named-suffix "~{named_suffix}" \
      --report-title "~{report_title}"

    ls -lh summary/results_by_phenotype/ || true
    wc -l summary/results_manifest.tsv summary/calibration_summary.tsv
  >>>

  output {
    File results_manifest = "summary/results_manifest.tsv"
    File calibration_summary = "summary/calibration_summary.tsv"
    File calibration_summary_md = "summary/calibration_summary.md"
    File lambda_gc_wide = "summary/lambda_gc_wide.tsv"
    File phewas_genomewide_hits = "summary/phewas_genomewide_hits.tsv"
    Array[File] results_tsvs_named = glob("summary/results_by_phenotype/*.felix.tsv")
    Array[File] qq_plots = glob("summary/qc/*/qq_joint_acpass.png")
    Array[File] manhattan_plots = glob("summary/qc/*/manhattan_joint.png")
    Array[File] top_hits_tables = glob("summary/qc/*/top_hits.tsv")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task CheckResumeInputs {
  input {
    Array[String] phenotypes
    Array[String] null_rda_names
    Array[String] variance_ratio_names
    Array[String] samples_used_names
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
  }

  command <<<
    set -euo pipefail
    PH=(~{sep=" " phenotypes})
    RDA=(~{sep=" " null_rda_names})
    VR=(~{sep=" " variance_ratio_names})
    SU=(~{sep=" " samples_used_names})
    n=${#PH[@]}
    if [[ "${n}" -lt 1 ]]; then echo "no phenotypes" >&2; exit 1; fi
    if [[ "${#RDA[@]}" -ne "${n}" || "${#VR[@]}" -ne "${n}" || "${#SU[@]}" -ne "${n}" ]]; then
      echo "array lengths differ: phenotypes=${n} rda=${#RDA[@]} varianceRatio=${#VR[@]} samples=${#SU[@]}" >&2
      exit 1
    fi
    for ((i = 0; i < n; i++)); do
      p="${PH[i]}"
      [[ "${RDA[i]}" == "${p}.null.rda" ]] || { echo "index ${i}: expected ${p}.null.rda, got ${RDA[i]}" >&2; exit 1; }
      [[ "${VR[i]}" == "${p}.varianceRatio.txt" ]] || { echo "index ${i}: expected ${p}.varianceRatio.txt, got ${VR[i]}" >&2; exit 1; }
      [[ "${SU[i]}" == "${p}.samples.txt" ]] || { echo "index ${i}: expected ${p}.samples.txt, got ${SU[i]}" >&2; exit 1; }
    done
    printf '%s\n' "${PH[@]}" > phenotypes.txt
    echo "OK: ${n} phenotypes line up with their null files"
  >>>

  output {
    Array[String] validated_phenotypes = read_lines("phenotypes.txt")
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}

workflow FelixPilotStep2Resume {
  input {
    # Outputs of the earlier run: call-Pack, call-MakeGRM, and call-Null/shard-N (one per phenotype).
    File packed_tar
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    Array[String] phenotypes
    Array[File] null_rdas
    Array[File] variance_ratios
    Array[File] samples_used

    # Same inputs as FelixPilot.
    File pheno_cov
    File run_felix_step2_script
    File summarize_script
    String chrom = "chr22"
    Int num_ancs = 5
    Int min_mac = 50
    Float pvalcutoff_of_haplotype = 0.05
    Float summarize_p_threshold = 0.00000005
    Int summarize_top_n = 50

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
  }

  scatter (f_rda in null_rdas) { String rda_name = basename(f_rda) }
  scatter (f_vr in variance_ratios) { String vr_name = basename(f_vr) }
  scatter (f_su in samples_used) { String su_name = basename(f_su) }

  call CheckResumeInputs as Check {
    input:
      phenotypes = phenotypes,
      null_rda_names = rda_name,
      variance_ratio_names = vr_name,
      samples_used_names = su_name,
      docker = docker
  }

  scatter (i in range(length(phenotypes))) {
    call RunFelixStep2 as Step2 {
      input:
        packed_tar = packed_tar,
        chrom = chrom,
        phenotype = Check.validated_phenotypes[i],
        null_rda = null_rdas[i],
        variance_ratio = variance_ratios[i],
        samples_used = samples_used[i],
        sparse_grm_mtx = sparse_grm_mtx,
        sparse_grm_sample_ids = sparse_grm_sample_ids,
        run_step2_script = run_felix_step2_script,
        num_ancs = num_ancs,
        min_mac = min_mac,
        pvalcutoff_of_haplotype = pvalcutoff_of_haplotype,
        docker = docker
    }
  }

  scatter (j in range(length(phenotypes))) {
    call ConcatPhenotypeScores as Concat {
      input:
        phenotype = Check.validated_phenotypes[j],
        shard_tsvs = [Step2.results_tsv[j]],
        docker = docker
    }
  }

  call SummarizeFelixResults as Summarize {
    input:
      result_tsvs = Concat.merged_tsv,
      phenotypes = Check.validated_phenotypes,
      pheno_cov = pheno_cov,
      summarize_script = summarize_script,
      p_threshold = summarize_p_threshold,
      top_n = summarize_top_n,
      docker = docker
  }

  output {
    Array[File] results_tsvs = Concat.merged_tsv
    Array[File] results_names = Concat.merged_name
    Array[File] results_tsvs_named = Summarize.results_tsvs_named
    File results_manifest = Summarize.results_manifest
    File calibration_summary = Summarize.calibration_summary
    File calibration_summary_md = Summarize.calibration_summary_md
    File lambda_gc_wide = Summarize.lambda_gc_wide
    File phewas_genomewide_hits = Summarize.phewas_genomewide_hits
    Array[File] qq_plots = Summarize.qq_plots
    Array[File] manhattan_plots = Summarize.manhattan_plots
    Array[File] top_hits_tables = Summarize.top_hits_tables
  }

  meta {
    description: "Temporary: FelixPilot Step2 -> concat -> summarize from the outputs of an earlier run."
    allowNestedInputs: true
  }
}
