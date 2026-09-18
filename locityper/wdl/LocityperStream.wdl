version 1.0

# Per-sample Locityper on a Nearline-efficient minicram.
# Upstream: EichlerLab/AoU_WDL locityper/locityper_stream.wdl (ValidateVariants).
# Extract: broadinstitute/str-analysis print_reads (CRAI → unique CRAM containers).
# Prep: notebooks/rw/locityper_00_prep_reference.ipynb
# Submit: notebooks/rw/locityper_01_run_stream.ipynb via `wb workflow`.

workflow LocityperStream {
    input {
        String sample_id

        File cram
        File crai

        File ref_fa_uncompressed
        File ref_fai_uncompressed

        File counts_jf
        File bed

        File locityper_db_tar_gz

        Int locityper_n_cpu = 32
        Int window_grab = 10000
        Int locityper_extra_mem_gb = 6
        Int n_preemptible = 2

        Int max_retry = 3
        Int wait_time = 30

        Int N_split = 200

        # Isaac's Terra run used illumina. PacBio HiFi: "hifi" / "pacbio" (see locityper --help).
        String technology = "illumina"

        # GRCh38 Locityper background interval; included in the minicram only.
        # Empty string skips the extra grab. CHM13 default is chr17:72950001-77450000.
        String bg_region_bed = "chr17\t72062001\t76562000"

        # Billing project for requester-pays GCS (AoU CRAMs). Empty → GOOGLE_CLOUD_PROJECT if set.
        String gcloud_project = ""

        Int minicram_disk_gb = 50
        Int minicram_mem_gb = 16
        Int minicram_n_preemptible = 0

        String locityper_docker = "eichlerlab/locityper:1.4.5.0"
        String print_reads_docker = "us-central1-docker.pkg.dev/broad-dsp-lrma/aou-lr/aou-locityper-print-reads:0.1.0"
        String util_docker = "python:3.11-slim"
    }

    call MakeMinicram {
        input:
            sample_id = sample_id,
            cram = cram,
            crai = crai,
            reference = ref_fa_uncompressed,
            bed = bed,
            window_grab = window_grab,
            bg_region_bed = bg_region_bed,
            gcloud_project = gcloud_project,
            max_retry = max_retry,
            wait_time = wait_time,
            disk_gb = minicram_disk_gb,
            mem_gb = minicram_mem_gb,
            n_preemptible = minicram_n_preemptible,
            docker = print_reads_docker
    }

    call LocityperPreprocess {
        input:
            sample_id = sample_id,
            minicram = MakeMinicram.minicram,
            minicrai = MakeMinicram.minicrai,
            reference = ref_fa_uncompressed,
            reference_index = ref_fai_uncompressed,
            counts_file = counts_jf,
            locityper_n_cpu = locityper_n_cpu,
            locityper_extra_mem_gb = locityper_extra_mem_gb,
            n_preemptible = n_preemptible,
            technology = technology,
            docker = locityper_docker
    }

    call SplitBed {
        input:
            bed = bed,
            N = N_split,
            docker = util_docker
    }

    scatter (split_bed in SplitBed.split_beds) {
        call LocityperGenotype {
            input:
                sample_id = sample_id,
                minicram = MakeMinicram.minicram,
                minicrai = MakeMinicram.minicrai,
                reference = ref_fa_uncompressed,
                reference_index = ref_fai_uncompressed,
                db_targz = locityper_db_tar_gz,
                preproc_tar_gz = LocityperPreprocess.preproc_tar_gz,
                bed = split_bed,
                locityper_n_cpu = locityper_n_cpu,
                locityper_extra_mem_gb = locityper_extra_mem_gb,
                n_preemptible = n_preemptible,
                docker = locityper_docker
        }
    }

    call CombineTarFiles {
        input:
            tar_files = LocityperGenotype.genotype_tar,
            sample_id = sample_id,
            docker = util_docker
    }

    call Summarize {
        input:
            sample_id = sample_id,
            genotype_tar = CombineTarFiles.combined_tar_gz,
            docker = util_docker
    }

    output {
        File summary_csv = Summarize.summary_csv
        File results_tar_gz = CombineTarFiles.combined_tar_gz
        File minicram = MakeMinicram.minicram
        File minicrai = MakeMinicram.minicrai
        File minicram_transfer_stats = MakeMinicram.transfer_stats
    }
}

task MakeMinicram {
    input {
        String sample_id
        File cram
        File crai
        File reference
        File bed
        Int window_grab
        String bg_region_bed
        String gcloud_project
        Int max_retry
        Int wait_time
        Int disk_gb
        Int mem_gb
        Int n_preemptible
        String docker
    }

    # Do not size disk from the WGS CRAM. localization_optional leaves it in Nearline;
    # Cromwell size() would still see the full object and over-allocate.
    parameter_meta {
        cram: { localization_optional: true }
    }

    String out_cram = sample_id + ".minicram.cram"

    command <<<
        set -euo pipefail
        date

        awk -v w="~{window_grab}" 'BEGIN { OFS="\t" } {
            s = $2 - w
            if (s < 0) s = 0
            print $1, s, $3 + w
        }' ~{bed} > intervals.bed

        if [ -n "~{bg_region_bed}" ]; then
            echo -e "~{bg_region_bed}" | awk -v w="~{window_grab}" 'BEGIN { OFS="\t" } {
                s = $2 - w
                if (s < 0) s = 0
                print $1, s, $3 + w
            }' >> intervals.bed
        fi

        echo "minicram intervals:"
        cat intervals.bed
        wc -l intervals.bed

        PROJ="~{gcloud_project}"
        if [ -z "${PROJ}" ]; then
            PROJ="${GOOGLE_CLOUD_PROJECT:-}"
        fi
        if [ -n "${PROJ}" ]; then
            echo "requester-pays project: ${PROJ}"
        fi

        run_print_reads() {
            local -a cmd=(
                python3 -m str_analysis.print_reads
                -R ~{reference}
                --read-index ~{crai}
                -L intervals.bed
                --padding 0
                -o ~{out_cram}
                --verbose
                --output-data-transfer-stats
            )
            if [ -n "${PROJ}" ]; then
                cmd+=(--gcloud-project "${PROJ}")
            fi
            cmd+=(~{cram})
            "${cmd[@]}"
        }

        MAX_RETRIES=~{max_retry}
        CURR_RETRIES=0
        WAIT_TIME=~{wait_time}
        until (( CURR_RETRIES == MAX_RETRIES )); do
            if run_print_reads; then
                break
            fi
            echo "print_reads failed; retry $(( CURR_RETRIES + 1 ))/${MAX_RETRIES}" >&2
            rm -f ~{out_cram} ~{out_cram}.crai *.data_transfer_stats.tsv
            CURR_RETRIES=$(( CURR_RETRIES + 1 ))
            if (( CURR_RETRIES == MAX_RETRIES )); then
                echo "print_reads failed after ${MAX_RETRIES} attempts" >&2
                exit 1
            fi
            sleep "${WAIT_TIME}"
        done

        if [ ! -s ~{out_cram} ]; then
            echo "print_reads wrote no CRAM" >&2
            exit 1
        fi
        if [ ! -s ~{out_cram}.crai ]; then
            echo "print_reads did not write ~{out_cram}.crai" >&2
            exit 1
        fi

        ls -lh ~{out_cram} ~{out_cram}.crai *.data_transfer_stats.tsv
        date
    >>>

    output {
        File minicram = out_cram
        File minicrai = out_cram + ".crai"
        File transfer_stats = glob("*.data_transfer_stats.tsv")[0]
    }

    runtime {
        memory: mem_gb + " GB"
        cpu: 4
        disks: "local-disk " + disk_gb + " HDD"
        preemptible: n_preemptible
        docker: docker
    }
}

task LocityperPreprocess {
    input {
        String sample_id
        File minicram
        File minicrai
        File reference
        File reference_index
        File counts_file
        Int locityper_n_cpu
        Int locityper_extra_mem_gb
        Int n_preemptible
        String technology
        String docker
    }

    Int disk_size = 4 + 3 * ceil(size([minicram, minicrai, counts_file, reference, reference_index], "GiB"))
    Int locityper_mem_gb = ceil(4.0 * locityper_n_cpu) + locityper_extra_mem_gb
    String preproc_tar = sample_id + ".preproc.tar.gz"

    command <<<
        set -euo pipefail
        date
        mv ~{reference} reference.fa
        mv ~{reference_index} reference.fa.fai
        mv ~{minicram} subset.cram
        mv ~{minicrai} subset.cram.crai

        locityper preproc -a subset.cram \
            -r reference.fa \
            -j ~{counts_file} \
            -@ ~{locityper_n_cpu} \
            --technology ~{technology} \
            -o locityper_preproc

        tar -czf ~{preproc_tar} locityper_preproc
        ls -lh ~{preproc_tar}
        date
    >>>

    output {
        File preproc_tar_gz = preproc_tar
    }

    runtime {
        memory: locityper_mem_gb + " GB"
        cpu: locityper_n_cpu
        disks: "local-disk " + disk_size + " HDD"
        preemptible: n_preemptible
        docker: docker
    }
}

task LocityperGenotype {
    input {
        File minicram
        File minicrai
        File reference
        File reference_index
        File db_targz
        File preproc_tar_gz
        String sample_id
        File bed
        Int locityper_n_cpu
        Int locityper_extra_mem_gb
        Int n_preemptible = 2
        String docker
    }

    Int disk_size = 4 + 3 * ceil(size([minicram, minicrai, reference, reference_index, db_targz, preproc_tar_gz, bed], "GiB"))
    Int locityper_mem_gb = ceil(4.0 * locityper_n_cpu) + locityper_extra_mem_gb
    String output_tar = sample_id + ".locityper.tar.gz"

    command <<<
        set -euo pipefail
        date
        mv ~{reference} reference.fa
        mv ~{reference_index} reference.fa.fai
        mv ~{minicram} subset.cram
        mv ~{minicrai} subset.cram.crai

        tar -xzf ~{preproc_tar_gz}
        tar -xzf ~{db_targz}

        mkdir -p out_dir/loci

        process_single_locus() {
            line="$1"
            locus_name=$(echo "$line" | cut -f4)
            mkdir -p "out_dir/loci/${locus_name}"

            locityper genotype -a subset.cram \
                -r reference.fa \
                -d vcf_db \
                -p locityper_preproc \
                --subset-loci "${locus_name}" \
                -o out_dir  >/dev/null 2>&1
        }
        export -f process_single_locus

        cat ~{bed} | /usr/bin/parallel --line-buffer -j ~{locityper_n_cpu} process_single_locus {}
        date
        find out_dir -type f -name "*.bam" -exec rm -f {} \;

        tar -czf ~{output_tar} out_dir
        date
    >>>

    output {
        File genotype_tar = output_tar
    }

    runtime {
        memory: locityper_mem_gb + " GB"
        cpu: locityper_n_cpu
        disks: "local-disk " + disk_size + " HDD"
        preemptible: n_preemptible
        docker: docker
    }
}

task SplitBed {
    input {
        File bed
        Int N = 20
        String docker
    }

    Int disk_size = 1 + ceil(size(bed, "GiB"))

    command <<<
        set -euo pipefail

        cat ~{bed} | split -l ~{N} - split_part_ && wc -l split_part_*
    >>>

    output {
        Array[File] split_beds = glob("split_part_*")
    }

    runtime {
        memory: "1 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: docker
    }
}

task CombineTarFiles {
    input {
        Array[File] tar_files
        String sample_id
        String docker
    }

    Int disk_size = 1 + 10*ceil(size(tar_files, "GiB"))
    String combined_tar = sample_id + ".combined.tar.gz"

    command <<<
        set -euo pipefail

        mkdir -p combined_temp

        for tar_file in ~{sep=" " tar_files}; do
            tar -xzf "$tar_file" -C combined_temp
        done

        tar -czf ~{combined_tar} -C combined_temp .

        rm -rf combined_temp
    >>>

    output {
        File combined_tar_gz = combined_tar
    }

    runtime {
        memory: "8 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: docker
    }
}

task Summarize {
    input {
        String sample_id
        File genotype_tar
        String docker
    }

    Int disk_size = 1 + 10*ceil(size(genotype_tar, "GiB"))

    command <<<
        set -euo pipefail

        tar -xzvf ~{genotype_tar}

        if [ -d out_dir ]; then
            mv out_dir ~{sample_id}
        elif [ ! -d ~{sample_id} ]; then
            echo "expected out_dir after extracting ~{genotype_tar}" >&2
            exit 1
        fi

        for dir in ~{sample_id}/loci/*/; do
            if [ ! -f "${dir}/res.json.gz" ]; then
                rm -rf "${dir}"
            fi
        done

        python3 - <<'PY'
import gzip
import json
import math
from pathlib import Path

sample_id = "~{sample_id}"
loci_root = Path(sample_id) / "loci"
out_path = Path("gts.filtered.csv")

rows = ["sample\tlocus\tgenotype\tquality\ttotal_reads\tunexpl_reads\tweight_dist\twarnings"]
if loci_root.is_dir():
    for entry in sorted(loci_root.iterdir(), key=lambda p: p.name):
        if not entry.is_dir():
            continue
        json_path = entry / "res.json.gz"
        if not json_path.is_file():
            continue
        with gzip.open(json_path, "rt") as fh:
            res = json.load(fh)
        line = f"{sample_id}\t{entry.name}\t"
        if "genotype" not in res:
            rows.append(line + "*")
            continue
        gt = res["genotype"]
        qual = math.floor(10 * float(res["quality"])) * 0.1
        total_reads = res.get("total_reads", "")
        unexpl_reads = res.get("unexpl_reads", "")
        weight_dist = res.get("weight_dist")
        weight_s = "" if weight_dist is None else f"{float(weight_dist):.5f}"
        warnings = ";".join(res.get("warnings", [])) or "*"
        rows.append(
            line + f"{gt}\t{qual:.1f}\t{total_reads}\t{unexpl_reads}\t{weight_s}\t{warnings}"
        )

out_path.write_text("\n".join(rows) + "\n")
print(f"wrote {len(rows) - 1} loci to {out_path}", flush=True)
PY
    >>>

    output {
        File summary_csv = "gts.filtered.csv"
    }

    runtime {
        memory: "4 GB"
        cpu: "1"
        disks: "local-disk " + disk_size + " HDD"
        preemptible: 1
        docker: docker
    }
}
