version 1.0

# Genome-wide FELIX for Terra Methods Repository (self-contained; no imports).
#
# Shape:
#   CheckChromVcfPairs (length-only; do not pass Array[File] to the check task)
#   BuildSaigePlinkAndSparseGRM once
#   FitFelixNull once per phenotype (fit_felix_null.R → null_export for RU scorer)
#   PackFelixla once per chromosome (joint_vcf OR phase_vcf + flare_vcf arrays)
#   RunFelixStep2 once per (phenotype × chromosome)
#   Concat → one merged TSV per phenotype
#   SummarizeFelixResults → phenotype-named copies, λGC / QQ / Manhattan
#
# Optional repeat-dosage branch:
#   ru_engine="felix" (default, needs repeat_catalog_bed):
#     AggregateRepeatLoci → WriteRuAdmixedVcf → RunFelixStep2RuVcf
#     One test per catalog repeat locus, REF-relative units, SPA + carrier QC.
#   ru_engine="rust" (needs simple_repeat_bed), comparison only:
#     AnnotateRuTest → ExtractRuDosage (dosage / collapse / split) → ScoreRuDosage
#   (tractor-mix-score --mode felix) per (phenotype × chromosome × encoding)
#
# Terra chrom-set launch: parallel Arrays
#   chroms      <- explicit chr names
#   flare_vcfs  <- FLARE anc VCFs (PackFelixla)
#   phase_vcfs  <- joint phased VCFs (PackFelixla + RU annotator)
# or joint_vcfs filling both phase and FLARE when GT:AN1:AN2 is in one file.

task CheckChromVcfPairs {
  input {
    Int n_chroms
    Int n_joint_vcfs
    Int n_phase_vcfs
    Int n_flare_vcfs
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
  }

  command <<<
    set -euo pipefail
    N_CHROM=~{n_chroms}
    N_JOINT=~{n_joint_vcfs}
    N_PHASE=~{n_phase_vcfs}
    N_FLARE=~{n_flare_vcfs}

    if [[ "${N_CHROM}" -lt 1 ]]; then
      echo "Need at least one chromosome" >&2
      exit 1
    fi

    if [[ "${N_JOINT}" -gt 0 ]]; then
      if [[ "${N_JOINT}" -ne "${N_CHROM}" ]]; then
        echo "joint_vcfs length (${N_JOINT}) != chroms length (${N_CHROM})" >&2
        exit 1
      fi
      echo "OK: ${N_CHROM} chrom × joint VCF pairs"
    else
      if [[ "${N_PHASE}" -ne "${N_CHROM}" || "${N_FLARE}" -ne "${N_CHROM}" ]]; then
        echo "phase_vcfs (${N_PHASE}) and flare_vcfs (${N_FLARE}) must both match chroms (${N_CHROM})" >&2
        exit 1
      fi
      echo "OK: ${N_CHROM} chrom × (phase + FLARE) pairs"
    fi
  >>>

  output {
    Int n_chroms_out = n_chroms
  }

  runtime {
    docker: docker
    cpu: 1
    memory: "1 GB"
    disks: "local-disk 10 HDD"
    preemptible: 3
  }
}

task BuildSaigePlinkAndSparseGRM {
  input {
    Array[File] grm_vcfs
    File analysis_samples
    File build_script
    File make_plink_keep_script
    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 300
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    chmod +x "~{build_script}"
    export MAKE_PLINK_KEEP_PY="~{make_plink_keep_script}"
    "~{build_script}" \
      --analysis-samples "~{analysis_samples}" \
      --relatedness-cutoff ~{relatedness_cutoff} \
      --n-threads ~{cpu} \
      --num-random-markers ~{num_random_markers} \
      --out-prefix saige_grm/plink \
      -- \
      ~{sep=" " grm_vcfs}

    cp -L saige_grm/sparseGRM.mtx saige_sparseGRM.mtx
    cp -L saige_grm/sparseGRM.sampleIDs.txt saige_sparseGRM.sampleIDs.txt
    cp saige_grm/plink.bed saige_plink.bed
    cp saige_grm/plink.bim saige_plink.bim
    cp saige_grm/plink.fam saige_plink.fam
    cp saige_grm/saige_sparse_grm.paths.tsv .
  >>>

  output {
    File plink_bed = "saige_plink.bed"
    File plink_bim = "saige_plink.bim"
    File plink_fam = "saige_plink.fam"
    File sparse_grm_mtx = "saige_sparseGRM.mtx"
    File sparse_grm_sample_ids = "saige_sparseGRM.sampleIDs.txt"
    File paths_tsv = "saige_sparse_grm.paths.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task FitFelixNull {
  input {
    File pheno_cov
    String phenotype
    File covariate_columns
    File analysis_samples
    File plink_bed
    File plink_bim
    File plink_fam
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    File fit_null_script
    String trait_type = "binary"
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 100
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    mkdir -p null
    ln -sf "~{plink_bed}" null/plink.bed
    ln -sf "~{plink_bim}" null/plink.bim
    ln -sf "~{plink_fam}" null/plink.fam

    Rscript "~{fit_null_script}" \
      --pheno-cov "~{pheno_cov}" \
      --phenotype "~{phenotype}" \
      --covariates "~{covariate_columns}" \
      --analysis-samples "~{analysis_samples}" \
      --plink-prefix null/plink \
      --sparse-grm "~{sparse_grm_mtx}" \
      --sparse-grm-ids "~{sparse_grm_sample_ids}" \
      --trait-type "~{trait_type}" \
      --step1-r /usr/local/bin/step1_fitNULLGLMM.R \
      --n-threads ~{cpu} \
      --out-prefix "null/~{phenotype}"

    cp "null/~{phenotype}.rda" "~{phenotype}.null.rda" || \
      cp $(ls -1 null/~{phenotype}*.rda | head -n1) "~{phenotype}.null.rda"
    cp "null/~{phenotype}.varianceRatio.txt" "~{phenotype}.varianceRatio.txt" || \
      cp $(ls -1 null/~{phenotype}*.varianceRatio.txt | head -n1) "~{phenotype}.varianceRatio.txt"
    cp "null/~{phenotype}.pheno.tsv" "~{phenotype}.pheno.tsv"
    cp "null/~{phenotype}.samples.txt" "~{phenotype}.samples.txt"
    cp "null/~{phenotype}.null_meta.tsv" "~{phenotype}.null_meta.tsv"

    NULL_EXPORT="null/~{phenotype}.null_export"
    if [[ ! -d "${NULL_EXPORT}" ]]; then
      echo "missing ${NULL_EXPORT} (fit_felix_null.R should export CSC layout)" >&2
      exit 1
    fi
    tar czf null_export.tar.gz -C null "$(basename "${NULL_EXPORT}")"
  >>>

  output {
    File null_rda = "~{phenotype}.null.rda"
    File variance_ratio = "~{phenotype}.varianceRatio.txt"
    File pheno_used = "~{phenotype}.pheno.tsv"
    File samples_used = "~{phenotype}.samples.txt"
    File null_meta = "~{phenotype}.null_meta.tsv"
    File null_export_tar = "null_export.tar.gz"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task PackFelixla {
  input {
    File phase_vcf
    File flare_vcf
    File analysis_samples
    String chrom
    Int num_ancs = 5
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 4
    Int memory_gb = 16
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 3.5
    Int preemptible = 0
  }

  Int disk_gb = ceil((size(phase_vcf, "GB") + size(flare_vcf, "GB")) * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p pack_vcf packed
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PREFIX="packed/${CHROM}"

    fix_vcf() {
      local in_vcf="$1"
      local out_vcf="$2"
      local hdr
      hdr=$(mktemp)
      bcftools view -h "${in_vcf}" \
        | awk '/^##FORMAT=<ID=GT,/{ if (gt++) next } { print }' \
        > "${hdr}"
      bcftools reheader -h "${hdr}" -o "${out_vcf}" "${in_vcf}"
      bcftools index -t "${out_vcf}" || bcftools index -c "${out_vcf}"
      rm -f "${hdr}"
    }

    if [[ "~{phase_vcf}" == "~{flare_vcf}" ]]; then
      fix_vcf "~{phase_vcf}" pack_vcf/joint.vcf.gz
      PHASE_IN=pack_vcf/joint.vcf.gz
      FLARE_IN=pack_vcf/joint.vcf.gz
    else
      fix_vcf "~{phase_vcf}" pack_vcf/phase.vcf.gz
      fix_vcf "~{flare_vcf}" pack_vcf/flare.vcf.gz
      PHASE_IN=pack_vcf/phase.vcf.gz
      FLARE_IN=pack_vcf/flare.vcf.gz
    fi

    felixla \
      --phase-vcf "${PHASE_IN}" \
      --flare-vcf "${FLARE_IN}" \
      --n-ancestries ~{num_ancs} \
      --keep "~{analysis_samples}" \
      --make-felixla \
      --out "${PREFIX}"

    if [[ ! -f "${PREFIX}.meta" || ! -f "${PREFIX}.samples" ]]; then
      echo "felixla did not write ${PREFIX}.meta / .samples" >&2
      ls -lh packed/ || true
      exit 1
    fi

    printf '%s\n' "${CHROM}" > packed/chrom.txt
    tar czf felixla_packed.tar.gz -C packed .
  >>>

  output {
    File packed_tar = "felixla_packed.tar.gz"
    String chrom_id = read_string("packed/chrom.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

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
    cat > shards.txt <<'EOF'
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
      --p-threshold ~{p_threshold}

    ls -lh summary/results_by_phenotype/ || true
    wc -l summary/results_manifest.tsv summary/calibration_summary.tsv
  >>>

  output {
    File results_manifest = "summary/results_manifest.tsv"
    File calibration_summary = "summary/calibration_summary.tsv"
    Array[File] results_tsvs_named = glob("summary/results_by_phenotype/*.felix.tsv")
    Array[File] qq_plots = glob("summary/qc/*/qq_cct.png")
    Array[File] manhattan_plots = glob("summary/qc/*/manhattan_cct.png")
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

task AnnotateRuTest {
  input {
    File source_vcf
    String chrom
    File simple_repeat_bed
    File annotate_script
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 2
    Int memory_gb = 8
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 2.5
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(source_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    OUT="annotated/annotated.ru.vcf.gz"
    mkdir -p annotated

    python3 "~{annotate_script}" \
      --vcf "~{source_vcf}" \
      --simple-repeat-bed "~{simple_repeat_bed}" \
      --out "${OUT}"

    bcftools index -t "${OUT}" || bcftools index -c "${OUT}"
    printf '%s\n' "${CHROM}" > annotated/chrom.txt
  >>>

  output {
    File annotated_vcf = "annotated/annotated.ru.vcf.gz"
    String chrom_id = read_string("annotated/chrom.txt")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task ExtractRuDosage {
  input {
    File annotated_vcf
    String chrom
    File analysis_samples
    Int num_ancs
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 4
    Int memory_gb = 16
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 3.5
    Int preemptible = 0
  }

  Int disk_gb = ceil(size(annotated_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    mkdir -p extract ordered
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")

    extract-tracts-flare \
      --vcf "~{annotated_vcf}" \
      --num-ancs ~{num_ancs} \
      --output-dir extract \
      --compress-output \
      --samples "~{analysis_samples}"

    python3 - <<'PY'
import re
import shutil
from pathlib import Path

src = Path("extract")
od = Path("ordered")
od.mkdir(parents=True, exist_ok=True)

for p in sorted(src.glob("*.dosage.txt.gz")):
    if ".collapse." in p.name or ".split." in p.name:
        continue
    m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
    if not m:
        raise SystemExit(f"unexpected dosage name: {p.name}")
    shutil.copy2(p, od / f"anc_{int(m.group(1)):02d}.dosage.txt.gz")
for p in sorted(src.glob("*.collapse.anc*.dosage.txt.gz")):
    m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
    if m:
        shutil.copy2(p, od / f"collapse_anc_{int(m.group(1)):02d}.dosage.txt.gz")
for p in sorted(src.glob("*.split.anc*.dosage.txt.gz")):
    m = re.search(r"\.anc(\d+)\.dosage\.txt\.gz$", p.name)
    if m:
        shutil.copy2(p, od / f"split_anc_{int(m.group(1)):02d}.dosage.txt.gz")

order = src / "dosage_sample_order.txt"
if not order.exists():
    raise SystemExit("missing extract/dosage_sample_order.txt")
shutil.copy2(order, od / "dosage_sample_order.txt")
PY

    cp "~{analysis_samples}" ordered/analysis_samples.txt
    python3 -c "import pathlib,sys; pathlib.Path('ordered/chrom.txt').write_text(sys.argv[1].strip().strip(chr(34)).strip(chr(39)) + chr(10))" "${CHROM}"
  >>>

  output {
    String chrom_id = read_string("ordered/chrom.txt")
    Array[File] dosage_files = glob("ordered/anc_*.dosage.txt.gz")
    Array[File] collapse_files = glob("ordered/collapse_anc_*.dosage.txt.gz")
    Array[File] split_files = glob("ordered/split_anc_*.dosage.txt.gz")
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

task ScoreRuDosage {
  input {
    String phenotype
    String chrom
    String encoding
    File null_export_tar
    Array[File] dosage_files
    Array[File] collapse_files
    Array[File] split_files
    Int score_threads = 8
    Int chunk_size = 2048
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 8
    Int memory_gb = 32
    Int disk_gb = 200
    Int preemptible = 0
  }

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    ENC=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{encoding}")
    OUT="${PHENO}.${CHROM}.ru_${ENC}.felix.tsv"

    tar xzf "~{null_export_tar}"
    NULL_DIR=$(tar tzf "~{null_export_tar}" | head -n1 | cut -d/ -f1)
    if [[ ! -d "${NULL_DIR}" ]]; then
      NULL_DIR=null_export
    fi

    DOSAGE_ARGS=()
    case "${ENC}" in
      dosage)
        for f in ~{sep=" " dosage_files}; do DOSAGE_ARGS+=("$f"); done
        ;;
      collapse)
        for f in ~{sep=" " collapse_files}; do DOSAGE_ARGS+=("$f"); done
        ;;
      split)
        for f in ~{sep=" " split_files}; do DOSAGE_ARGS+=("$f"); done
        ;;
      *)
        echo "encoding must be dosage, collapse, or split; got ${ENC}" >&2
        exit 1
        ;;
    esac

    if [[ "${#DOSAGE_ARGS[@]}" -lt 1 ]]; then
      echo "no dosage files for encoding=${ENC}" >&2
      exit 1
    fi

    tractor-mix-score \
      --null-export "${NULL_DIR}" \
      --dosage-files "${DOSAGE_ARGS[@]}" \
      --out "${OUT}" \
      --mode felix \
      --threads ~{score_threads} \
      --chunk-size ~{chunk_size}

    cp "${OUT}" results.ru.felix.tsv
    printf '%s\n' "${OUT}" > results.output_name.txt
  >>>

  output {
    File results_tsv = "results.ru.felix.tsv"
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

# One tested unit per catalog repeat locus: every length-changing record inside a
# TRExplorer interval contributes its signed length change to that haplotype's
# locus dosage (scripts/aggregate_repeat_loci.py --out-vcf). This replaces the
# record-level annotator for ru_engine="felix" — see felix/REPEAT_DOSAGE.md §6.
task AggregateRepeatLoci {
  input {
    File source_vcf
    File? source_vcf_index
    String chrom
    File repeat_catalog_bed
    File analysis_samples
    File aggregate_script
    Int num_ancs = 5
    Int flank = 10
    Int max_bp = 100000
    String apply_filters = "PASS,."
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 4
    Int memory_gb = 16
    Int disk_gb_floor = 100
    Float disk_gb_multiplier = 3.0
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(source_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    mkdir -p loci src

    # The aggregator seeks with bcftools, so the source VCF needs an index next
    # to it. Cromwell localizes the index (when given) to a different directory,
    # so symlink both into one place; index here when none was provided.
    SRC="src/$(basename "~{source_vcf}")"
    ln -s "~{source_vcf}" "${SRC}"
    IDX="~{default='' source_vcf_index}"
    if [[ -n "${IDX}" ]]; then
      ln -s "${IDX}" "src/$(basename "${IDX}")"
    else
      bcftools index -t "${SRC}" || bcftools index -c "${SRC}"
    fi

    # --num-ancs fails loudly if FLARE emitted a code the downstream tests
    # cannot represent, instead of silently folding it into another ancestry.
    python3 "~{aggregate_script}" \
      --vcf "${SRC}" \
      --catalog-bed "~{repeat_catalog_bed}" \
      --chrom "${CHROM}" \
      --samples "~{analysis_samples}" \
      --apply-filters "~{apply_filters}" \
      --flank ~{flank} \
      --max-bp ~{max_bp} \
      --with-ancestry \
      --num-ancs ~{num_ancs} \
      --out-vcf loci/loci.vcf \
      --out-alleles loci/loci.alleles.tsv.gz \
      --out-summary loci/loci.summary.json

    bcftools view -Oz -o loci/loci.vcf.gz loci/loci.vcf
    bcftools index -c loci/loci.vcf.gz
    rm -f loci/loci.vcf
    cat loci/loci.summary.json
  >>>

  output {
    File locus_vcf = "loci/loci.vcf.gz"
    File locus_vcf_index = "loci/loci.vcf.gz.csi"
    File alleles_tsv = "loci/loci.alleles.tsv.gz"
    File summary_json = "loci/loci.summary.json"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

# Repeat-SV dosages for FELIX's admixed VCF path: DS{k} / ANC{k} from RU_TEST
# loci (REF-relative repeat units, scaled per locus; scripts/write_admixed_dosage_vcf.py).
task WriteRuAdmixedVcf {
  input {
    File annotated_vcf
    String chrom
    File analysis_samples
    File write_script
    Int num_ancs = 5
    String missing_policy = "ref"
    Float max_hap_missing = 0.1
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 2
    Int memory_gb = 8
    Int disk_gb_floor = 50
    Float disk_gb_multiplier = 3.0
    Int preemptible = 1
  }

  Int disk_gb = ceil(size(annotated_vcf, "GB") * disk_gb_multiplier) + disk_gb_floor

  command <<<
    set -euo pipefail
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    mkdir -p ru
    python3 "~{write_script}" \
      --vcf "~{annotated_vcf}" \
      --num-ancs ~{num_ancs} \
      --samples "~{analysis_samples}" \
      --missing "~{missing_policy}" \
      --max-hap-missing ~{max_hap_missing} \
      --out ru/ru.admixed.vcf \
      --stats-json ru/ru.admixed.stats.json
    bcftools view -Oz -o ru/ru.admixed.vcf.gz ru/ru.admixed.vcf
    bcftools index -c ru/ru.admixed.vcf.gz
    rm -f ru/ru.admixed.vcf
    cat ru/ru.admixed.stats.json
    printf '%s\n' "${CHROM}" > ru/chrom.txt
  >>>

  output {
    File admixed_vcf = "ru/ru.admixed.vcf.gz"
    File admixed_vcf_index = "ru/ru.admixed.vcf.gz.csi"
    File stats_json = "ru/ru.admixed.stats.json"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

# FELIX step 2 on the repeat dosage VCF with carrier-count QC
# (FELIX_DOSAGE_QC=carrier; needs the patched image, felix/patches/).
task RunFelixStep2RuVcf {
  input {
    File admixed_vcf
    File admixed_vcf_index
    String chrom
    String phenotype
    File null_rda
    File variance_ratio
    File samples_used
    File sparse_grm_mtx
    File sparse_grm_sample_ids
    File run_step2_script
    Int num_ancs = 5
    Int min_carriers = 20
    Float pvalcutoff_of_haplotype = 0.05
    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
    Int cpu = 2
    Int memory_gb = 16
    Int disk_gb = 50
    Int preemptible = 1
  }

  command <<<
    set -euo pipefail
    mkdir -p step2
    CHROM=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{chrom}")
    PHENO=$(python3 -c 'import sys; print(sys.argv[1].strip().strip(chr(34)).strip(chr(39)))' "~{phenotype}")
    ln -sf "~{null_rda}" "step2/${PHENO}.rda"
    ln -sf "~{variance_ratio}" "step2/${PHENO}.varianceRatio.txt"
    ln -sf "~{admixed_vcf}" ru.admixed.vcf.gz
    ln -sf "~{admixed_vcf_index}" ru.admixed.vcf.gz.csi

    # --min-mac is a carrier count here (samples with a non-zero dosage).
    Rscript "~{run_step2_script}" \
      --vcf-file ru.admixed.vcf.gz \
      --vcf-index ru.admixed.vcf.gz.csi \
      --dosage-qc carrier \
      --chrom "${CHROM}" \
      --null-prefix "step2/${PHENO}" \
      --sample-file "~{samples_used}" \
      --sparse-grm "~{sparse_grm_mtx}" \
      --sparse-grm-ids "~{sparse_grm_sample_ids}" \
      --min-mac ~{min_carriers} \
      --n-ancestries ~{num_ancs} \
      --pvalcutoff-of-haplotype ~{pvalcutoff_of_haplotype} \
      --n-threads 1 \
      --out-tsv "${PHENO}.ru.felix.tsv" \
      --out-raw "${PHENO}.ru.felix.raw.txt"

    cp "${PHENO}.ru.felix.tsv" results.ru.felix.tsv
  >>>

  output {
    File results_tsv = "results.ru.felix.tsv"
  }

  runtime {
    docker: docker
    cpu: cpu
    memory: memory_gb + " GB"
    disks: "local-disk " + disk_gb + " HDD"
    preemptible: preemptible
  }
}

workflow FelixGenome {
  input {
    Array[String] chroms

  # Parallel per-chrom VCFs (same length / order as chroms).
  # Supply joint_vcfs OR both phase_vcfs and flare_vcfs.
    Array[File]? joint_vcfs
    Array[File]? phase_vcfs
    Array[File]? flare_vcfs
    Int num_ancs = 5

    File analysis_samples
    File pheno_cov
    File selected_phenotypes
    File covariate_columns
    String trait_type = "binary"

    Array[File] grm_vcfs

    File build_saige_grm_script
    File make_plink_keep_script
    File fit_felix_null_script
    File run_felix_step2_script
    File summarize_script

  # Optional repeat-dosage branch. Both engines read GT:AN1:AN2, so pass
  # joint_vcfs (PropagateFlareAncestry output), not phase_vcfs alone.
  #
  # "felix" (default): one test per catalog repeat locus. Needs
  #   repeat_catalog_bed (TRExplorer intervals: chrom start end TRID motifs),
  #   aggregate_repeat_loci_script and write_admixed_dosage_script, plus the
  #   patched image for carrier QC. felix/REPEAT_DOSAGE.md.
  # "rust": legacy record-level tractor-mix-score encodings (no SPA, one test
  #   per VCF record), for the encoding comparison only. Needs
  #   simple_repeat_bed and annotate_repeat_units_script.
    String ru_engine = "felix"

    File? repeat_catalog_bed
    File? aggregate_repeat_loci_script
    File? write_admixed_dosage_script
    Int ru_min_carriers = 20
    Int ru_flank = 10
    Int ru_max_bp = 100000
    String ru_apply_filters = "PASS,."
    String ru_missing_policy = "ref"
    Float ru_max_hap_missing = 0.1

    File? simple_repeat_bed
    File? annotate_repeat_units_script

    Float relatedness_cutoff = 0.05
    Int num_random_markers = 2000
    Int min_mac = 50
    Float pvalcutoff_of_haplotype = 0.05
    Float summarize_p_threshold = 0.00000005
    Int summarize_top_n = 50

    Int pack_cpu = 4
    Int pack_memory_gb = 16
    Int pack_disk_gb_floor = 100
    Float pack_disk_gb_multiplier = 3.5
    Int step2_disk_gb_floor = 50
    Float step2_disk_gb_multiplier = 4.0
    Int ru_disk_gb_floor = 100
    Float ru_disk_gb_multiplier = 3.5
    Int score_threads = 8
    Int chunk_size = 2048
    Int score_disk_gb = 200

    String docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/felix-pilot:0.2.0"
  }

  Array[String] phenotypes = read_lines(selected_phenotypes)
  Array[String] ru_encodings = ["dosage", "collapse", "split"]

  call CheckChromVcfPairs as Check {
    input:
      n_chroms = length(chroms),
      n_joint_vcfs = length(select_first([joint_vcfs, []])),
      n_phase_vcfs = length(select_first([phase_vcfs, []])),
      n_flare_vcfs = length(select_first([flare_vcfs, []])),
      docker = docker
  }

  call BuildSaigePlinkAndSparseGRM as MakeGRM {
    input:
      grm_vcfs = grm_vcfs,
      analysis_samples = analysis_samples,
      build_script = build_saige_grm_script,
      make_plink_keep_script = make_plink_keep_script,
      relatedness_cutoff = relatedness_cutoff,
      num_random_markers = num_random_markers,
      docker = docker
  }

  scatter (pheno in phenotypes) {
    call FitFelixNull as Null {
      input:
        pheno_cov = pheno_cov,
        phenotype = pheno,
        covariate_columns = covariate_columns,
        analysis_samples = analysis_samples,
        plink_bed = MakeGRM.plink_bed,
        plink_bim = MakeGRM.plink_bim,
        plink_fam = MakeGRM.plink_fam,
        sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
        sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
        fit_null_script = fit_felix_null_script,
        trait_type = trait_type,
        docker = docker
    }
  }

  scatter (i in range(length(chroms))) {
    String chrom = chroms[i]
    # Cromwell cannot index an optional array; select_first unwraps it
    # (CheckChromVcfPairs has already checked which inputs were given).
    File phase_for_pack = if defined(joint_vcfs) then select_first([joint_vcfs])[i] else select_first([phase_vcfs])[i]
    File flare_for_pack = if defined(joint_vcfs) then select_first([joint_vcfs])[i] else select_first([flare_vcfs])[i]
    File phase_for_ru = if defined(joint_vcfs) then select_first([joint_vcfs])[i] else select_first([phase_vcfs])[i]

    call PackFelixla as Pack {
      input:
        phase_vcf = phase_for_pack,
        flare_vcf = flare_for_pack,
        analysis_samples = analysis_samples,
        chrom = chrom,
        num_ancs = num_ancs,
        docker = docker,
        cpu = pack_cpu,
        memory_gb = pack_memory_gb,
        disk_gb_floor = pack_disk_gb_floor,
        disk_gb_multiplier = pack_disk_gb_multiplier
    }

    # Locus-level dosage (default). One test per catalog repeat locus: the
    # record-level annotator is not used here, because it tests each VCF record
    # separately and drops records whose alleles are not a whole number of
    # repeat units (felix/REPEAT_DOSAGE.md §6.2, felix/eval/record_vs_locus/).
    if (ru_engine == "felix" && defined(repeat_catalog_bed)
        && defined(aggregate_repeat_loci_script) && defined(write_admixed_dosage_script)) {
      call AggregateRepeatLoci as RuLoci {
        input:
          source_vcf = phase_for_ru,
          chrom = chrom,
          repeat_catalog_bed = select_first([repeat_catalog_bed]),
          analysis_samples = analysis_samples,
          aggregate_script = select_first([aggregate_repeat_loci_script]),
          num_ancs = num_ancs,
          flank = ru_flank,
          max_bp = ru_max_bp,
          apply_filters = ru_apply_filters,
          docker = docker,
          disk_gb_floor = ru_disk_gb_floor,
          disk_gb_multiplier = ru_disk_gb_multiplier
      }

      call WriteRuAdmixedVcf as RuVcf {
        input:
          annotated_vcf = RuLoci.locus_vcf,
          chrom = chrom,
          analysis_samples = analysis_samples,
          write_script = select_first([write_admixed_dosage_script]),
          num_ancs = num_ancs,
          missing_policy = ru_missing_policy,
          max_hap_missing = ru_max_hap_missing,
          docker = docker
      }

      scatter (j in range(length(phenotypes))) {
        call RunFelixStep2RuVcf as RuStep2 {
          input:
            admixed_vcf = RuVcf.admixed_vcf,
            admixed_vcf_index = RuVcf.admixed_vcf_index,
            chrom = chrom,
            phenotype = phenotypes[j],
            null_rda = Null.null_rda[j],
            variance_ratio = Null.variance_ratio[j],
            samples_used = Null.samples_used[j],
            sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
            sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
            run_step2_script = run_felix_step2_script,
            num_ancs = num_ancs,
            min_carriers = ru_min_carriers,
            pvalcutoff_of_haplotype = pvalcutoff_of_haplotype,
            docker = docker
        }
      }
    }

    # Legacy record-level encodings, kept for the encoding comparison only.
    if (ru_engine == "rust" && defined(simple_repeat_bed)
        && defined(annotate_repeat_units_script)) {
      call AnnotateRuTest as RuAnnot {
        input:
          source_vcf = phase_for_ru,
          chrom = chrom,
          simple_repeat_bed = select_first([simple_repeat_bed]),
          annotate_script = select_first([annotate_repeat_units_script]),
          docker = docker
      }

      call ExtractRuDosage as RuExtract {
        input:
          annotated_vcf = RuAnnot.annotated_vcf,
          chrom = chrom,
          analysis_samples = analysis_samples,
          num_ancs = num_ancs,
          docker = docker,
          disk_gb_floor = ru_disk_gb_floor,
          disk_gb_multiplier = ru_disk_gb_multiplier
      }

      scatter (j in range(length(phenotypes))) {
        scatter (enc in ru_encodings) {
          call ScoreRuDosage as RuScore {
            input:
              phenotype = phenotypes[j],
              chrom = chrom,
              encoding = enc,
              null_export_tar = Null.null_export_tar[j],
              dosage_files = RuExtract.dosage_files,
              collapse_files = RuExtract.collapse_files,
              split_files = RuExtract.split_files,
              score_threads = score_threads,
              chunk_size = chunk_size,
              docker = docker,
              disk_gb = score_disk_gb
          }
        }
      }
    }

    scatter (j in range(length(phenotypes))) {
      call RunFelixStep2 as Step2 {
        input:
          packed_tar = Pack.packed_tar,
          chrom = chrom,
          phenotype = phenotypes[j],
          null_rda = Null.null_rda[j],
          variance_ratio = Null.variance_ratio[j],
          samples_used = Null.samples_used[j],
          sparse_grm_mtx = MakeGRM.sparse_grm_mtx,
          sparse_grm_sample_ids = MakeGRM.sparse_grm_sample_ids,
          run_step2_script = run_felix_step2_script,
          num_ancs = num_ancs,
          min_mac = min_mac,
          pvalcutoff_of_haplotype = pvalcutoff_of_haplotype,
          docker = docker,
          disk_gb_floor = step2_disk_gb_floor,
          disk_gb_multiplier = step2_disk_gb_multiplier
      }
    }
  }

  Array[Array[File]] results_by_pheno = transpose(Step2.results_tsv)
  # Chroms with RU results (felix engine), then per phenotype across chroms.
  Array[Array[File]] ru_by_chrom = select_all(RuStep2.results_tsv)
  Boolean has_ru = length(ru_by_chrom) > 0

  if (has_ru) {
    Array[Array[File]] ru_by_pheno = transpose(ru_by_chrom)
    scatter (i in range(length(phenotypes))) {
      call ConcatPhenotypeScores as RuConcat {
        input:
          phenotype = phenotypes[i] + ".ru",
          shard_tsvs = ru_by_pheno[i],
          docker = docker
      }
    }
  }

  scatter (i in range(length(phenotypes))) {
    call ConcatPhenotypeScores as Concat {
      input:
        phenotype = phenotypes[i],
        shard_tsvs = results_by_pheno[i],
        docker = docker
    }
  }

  call SummarizeFelixResults as Summarize {
    input:
      result_tsvs = Concat.merged_tsv,
      phenotypes = phenotypes,
      pheno_cov = pheno_cov,
      summarize_script = summarize_script,
      p_threshold = summarize_p_threshold,
      top_n = summarize_top_n,
      docker = docker
  }

  output {
    Int n_chroms = Check.n_chroms_out
    File plink_bed = MakeGRM.plink_bed
    File sparse_grm_mtx = MakeGRM.sparse_grm_mtx
    Array[File] null_rdas = Null.null_rda
    Array[File] null_export_tars = Null.null_export_tar
    Array[File] results_tsvs = Concat.merged_tsv
    Array[File] results_names = Concat.merged_name
    Array[File] results_tsvs_named = Summarize.results_tsvs_named
    File results_manifest = Summarize.results_manifest
    File calibration_summary = Summarize.calibration_summary
    Array[File] qq_plots = Summarize.qq_plots
    Array[File] manhattan_plots = Summarize.manhattan_plots
    Array[File] top_hits_tables = Summarize.top_hits_tables
    Array[Array[File]] results_tsvs_by_chrom = Step2.results_tsv
    Array[File] packed_tars = Pack.packed_tar
    Array[String] chrom_ids = Pack.chrom_id
    Array[File]? ru_results_tsvs = RuConcat.merged_tsv
    Array[File?] ru_admixed_vcfs = RuVcf.admixed_vcf
    Array[File?] ru_admixed_stats = RuVcf.stats_json
    # Locus/haplotype tallies travel with the results: missingness, possible
    # double counts between the indel and SV calls, ancestry inconsistency.
    Array[File?] ru_locus_summaries = RuLoci.summary_json
    Array[File?] ru_locus_alleles = RuLoci.alleles_tsv
    Array[File?] ru_locus_vcfs = RuLoci.locus_vcf
    # ru_engine="rust" only: per-chrom, per-phenotype, per-encoding score TSVs
    # for the encoding comparison. Not merged into Summarize.
    Array[Array[Array[File]]?] ru_rust_encoding_tsvs = RuScore.results_tsv
  }

  meta {
    description: "Genome-wide FELIX: shared SAIGE mtx GRM + FELIX nulls; per-chr felixla pack/step2; optional locus-level repeat dosages through FELIX step 2 (carrier QC)."
    allowNestedInputs: true
  }
}
