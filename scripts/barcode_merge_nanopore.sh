#!/usr/bin/env bash
# Merge per-barcode Nanopore BAMs into one sorted, indexed BAM per barcode.
#
# Activate your samtools env first (e.g. `conda activate python_collection`),
# then run interactively or submit with sbatch.
#
# Interactive:   ./barcode_merge_nanopore.sh -i /path/to/bam_pass -t 8
# SLURM:         sbatch barcode_merge_nanopore.sh   (edit INPUT_DIR below first)

#SBATCH --job-name=merge_bams
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=0
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

# When submitted with sbatch, this is the bam_pass dir we work on.
INPUT_DIR="/path/to/bam_pass"

# CLI flags override the defaults above.
THREADS=""
OUTPUT_DIR=""
FORCE=0

while getopts "i:o:t:fh" opt; do
    case "$opt" in
        i) INPUT_DIR="$OPTARG" ;;
        o) OUTPUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        f) FORCE=1 ;;
        h) sed -n '2,9p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) exit 1 ;;
    esac
done

# Threads: take from SLURM if we're inside a job, else default to 4.
[[ -z "$THREADS" ]] && THREADS="${SLURM_CPUS_PER_TASK:-4}"

shift $((OPTIND - 1))
[[ $# -gt 0 && -d "$1" ]] && INPUT_DIR="$1"

INPUT_DIR="$(cd "$INPUT_DIR" && pwd)"
[[ -z "$OUTPUT_DIR" ]] && OUTPUT_DIR="$INPUT_DIR/merged_bam"
mkdir -p "$OUTPUT_DIR"

if ! command -v samtools >/dev/null 2>&1; then
    echo "samtools not found on PATH — activate the right conda env first." >&2
    exit 127
fi

# Banner so you know exactly what this run is doing.
echo "=================================================================="
echo " Barcode merge — Nanopore BAMs"
echo "=================================================================="
echo " host         : $(hostname)"
echo " slurm job    : ${SLURM_JOB_ID:-(not in slurm)}"
echo " started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo " samtools     : $(samtools --version | head -n1)"
echo " input  dir   : $INPUT_DIR"
echo " output dir   : $OUTPUT_DIR"
echo " threads      : $THREADS"
echo " force        : $FORCE"
echo "------------------------------------------------------------------"
echo

# Collect the folders we want to process — barcode01..barcode99 and unclassified.
shopt -s nullglob
folders=()
for d in "$INPUT_DIR"/barcode[0-9][0-9] "$INPUT_DIR"/unclassified; do
    [[ -d "$d" ]] && folders+=("$d")
done
shopt -u nullglob

if [[ ${#folders[@]} -eq 0 ]]; then
    echo "Nothing to do — no barcodeNN/ or unclassified/ folders under $INPUT_DIR" >&2
    exit 1
fi

# Quick pre-flight summary so we know what's ahead.
echo "Found ${#folders[@]} folder(s) to process:"
total_input=0
for d in "${folders[@]}"; do
    n=$(find "$d" -maxdepth 1 -type f -name '*.bam' | wc -l)
    sz=$(du -sh "$d" 2>/dev/null | cut -f1)
    printf "  %-15s  %3d .bam files   %s\n" "$(basename "$d")" "$n" "$sz"
    total_input=$((total_input + n))
done
echo "  ---------------------------------------"
printf "  %-15s  %3d .bam files total\n" "TOTAL" "$total_input"
echo

# Process each folder. Pipe merge -> sort so no big intermediate file hits disk.
job_start=$SECONDS
fail=0
done_count=0
skip_count=0

for d in "${folders[@]}"; do
    name="$(basename "$d")"
    out_bam="$OUTPUT_DIR/${name}.bam"
    out_bai="${out_bam}.bai"

    mapfile -t bams < <(find "$d" -maxdepth 1 -type f -name '*.bam' | sort)

    if [[ ${#bams[@]} -eq 0 ]]; then
        echo ">> $name: no .bam files, skipping"
        skip_count=$((skip_count + 1))
        echo
        continue
    fi

    if [[ -s "$out_bam" && -s "$out_bai" && "$FORCE" -ne 1 ]]; then
        echo ">> $name: already merged ($(du -h "$out_bam" | cut -f1)), use -f to overwrite — skipping"
        skip_count=$((skip_count + 1))
        echo
        continue
    fi

    step_start=$SECONDS
    echo ">> $name"
    echo "   merging ${#bams[@]} bam(s)  -->  $out_bam"

    if ! samtools merge -@ "$THREADS" -u -f -O bam - "${bams[@]}" \
            | samtools sort -@ "$THREADS" -o "$out_bam" - ; then
        echo "   !! merge|sort FAILED for $name" >&2
        fail=$((fail + 1))
        echo
        continue
    fi

    echo "   indexing..."
    if ! samtools index -@ "$THREADS" "$out_bam"; then
        echo "   !! index FAILED for $name" >&2
        fail=$((fail + 1))
        echo
        continue
    fi

    # quick stats for the merged file
    reads=$(samtools view -@ "$THREADS" -c "$out_bam")
    size=$(du -h "$out_bam" | cut -f1)
    elapsed=$((SECONDS - step_start))
    printf "   done — %s reads, %s, %dm%ds\n" "$reads" "$size" $((elapsed / 60)) $((elapsed % 60))
    done_count=$((done_count + 1))
    echo
done

# Final summary.
total_elapsed=$((SECONDS - job_start))
echo "=================================================================="
echo " Summary"
echo "------------------------------------------------------------------"
echo " merged       : $done_count folder(s)"
echo " skipped      : $skip_count folder(s)"
echo " failed       : $fail folder(s)"
echo " total time   : $((total_elapsed / 60))m $((total_elapsed % 60))s"
echo " output dir   : $OUTPUT_DIR"
echo " finished     : $(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================================="

[[ "$fail" -ne 0 ]] && exit 1 || exit 0
