#!/usr/bin/env bash
# Filter each merged barcode BAM down to primary mapped reads only (-F 2308),
# then index the result. The output BAMs are what you want for everything
# downstream of DCS counting: bigwigs, read-start BEDs, IGV, anything that
# shouldn't double-count split/multi-aligning reads.
#
# Run AFTER count_dcs_spikein.sh (which needs the unmapped reads for DCS
# detection). The original merged BAMs are untouched — this writes a new
# folder of filtered BAMs alongside.
#
# `samtools view -F 2308` drops:
#     4    = unmapped
#     256  = secondary alignment
#     2048 = supplementary alignment
# i.e. you keep exactly one primary alignment per mapped read.
#
# Filtering preserves sort order, so we just view + index — no need to re-sort.
#
# Activate your env first so samtools is on PATH, then:
#     ./filter_primary_bams.sh -b /path/to/merged_bam -o /path/to/primary_bams -t 8
#  or sbatch filter_primary_bams.sh   (edit BAM_DIR below first)

#SBATCH --job-name=filter_primary
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=0
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

# Default — edit for sbatch use, or override with -b
BAM_DIR="/path/to/merged_bam"
OUT_DIR=""
THREADS=""
FORCE=0
INCLUDE_UNCLASSIFIED=0

while getopts "b:o:t:fuh" opt; do
    case "$opt" in
        b) BAM_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        f) FORCE=1 ;;
        u) INCLUDE_UNCLASSIFIED=1 ;;
        h) sed -n '2,22p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) exit 1 ;;
    esac
done

shift $((OPTIND - 1))
[[ $# -gt 0 && -d "$1" ]] && BAM_DIR="$1"

[[ -z "$THREADS" ]] && THREADS="${SLURM_CPUS_PER_TASK:-4}"

# Validate numeric thread value
if ! [[ "$THREADS" =~ ^[0-9]+$ ]] || [[ "$THREADS" -le 0 ]]; then
    echo "Error: THREADS must be a positive integer, got: $THREADS" >&2
    exit 1
fi

BAM_DIR="$(cd "$BAM_DIR" && pwd)" || { echo "Error: BAM_DIR does not exist: $BAM_DIR" >&2; exit 1; }
[[ -z "$OUT_DIR" ]] && OUT_DIR="$(dirname "$BAM_DIR")/primary_bams"
mkdir -p "$OUT_DIR"

if ! command -v samtools >/dev/null 2>&1; then
    echo "samtools not on PATH — activate your env." >&2
    exit 127
fi

echo "=================================================================="
echo " filter to primary mapped reads (-F 2308)"
echo "=================================================================="
echo " host         : $(hostname)"
echo " slurm job    : ${SLURM_JOB_ID:-(not in slurm)}"
echo " started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo " samtools     : $(samtools --version | head -n1)"
echo " input dir    : $BAM_DIR"
echo " output dir   : $OUT_DIR"
echo " threads      : $THREADS"
echo " force        : $FORCE"
echo " unclassified : $([[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && echo "included" || echo "skipped")"
echo "------------------------------------------------------------------"
echo

# IMPROVED: Flexible barcode pattern matching (supports any number of digits)
shopt -s nullglob extglob
candidates=("$BAM_DIR"/barcode+([0-9]).bam)
[[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && candidates+=("$BAM_DIR"/unclassified.bam)
shopt -u nullglob extglob
bams=()
for b in "${candidates[@]}"; do [[ -f "$b" ]] && bams+=("$b"); done

if [[ ${#bams[@]} -eq 0 ]]; then
    echo "No barcode BAMs in $BAM_DIR" >&2
    exit 1
fi

echo "Found ${#bams[@]} BAM(s) to filter."
echo

job_start=$SECONDS
done_count=0
skip_count=0
fail=0

for b in "${bams[@]}"; do
    name="$(basename "$b" .bam)"
    out="$OUT_DIR/${name}.bam"
    bai="${out}.bai"

    if [[ -s "$out" && -s "$bai" && "$FORCE" -ne 1 ]]; then
        echo ">> $name: already filtered ($(du -h "$out" | cut -f1)) — skipping (use -f to overwrite)"
        skip_count=$((skip_count + 1))
        echo
        continue
    fi

    step_start=$SECONDS
    echo ">> $name"
    in_reads=$(samtools view -@ "$THREADS" -c "$b")
    echo "   input reads   : $in_reads"
    echo -n "   filtering -F 2308 -> $out..."

    if ! samtools view -@ "$THREADS" -b -F 2308 -o "$out" "$b" 2>/dev/null; then
        echo " FAILED"
        fail=$((fail + 1))
        echo
        continue
    fi
    if ! samtools index -@ "$THREADS" "$out" 2>/dev/null; then
        echo " indexed FAILED"
        fail=$((fail + 1))
        echo
        continue
    fi

    out_reads=$(samtools view -@ "$THREADS" -c "$out")
    dropped=$((in_reads - out_reads))
    elapsed=$((SECONDS - step_start))
    size=$(du -h "$out" | cut -f1)
    # IMPROVED: Safer arithmetic with validation
    if [[ $in_reads -gt 0 ]]; then
        pct=$(awk -v d="$dropped" -v t="$in_reads" 'BEGIN{ printf "%.2f", 100*d/t}')
    else
        pct="0"
    fi
    printf " done\n   kept %s reads, dropped %s (%s%%), %s, %dm %ds\n\n" \
        "$out_reads" "$dropped" "$pct" "$size" $((elapsed/60)) $((elapsed%60))
    done_count=$((done_count + 1))
done

total=$((SECONDS - job_start))
echo "=================================================================="
echo " Summary"
echo "------------------------------------------------------------------"
printf " filtered   : %d\n" "$done_count"
printf " skipped    : %d\n" "$skip_count"
printf " failed     : %d\n" "$fail"
printf " total time : %dm %ds\n" $((total/60)) $((total%60))
printf " output dir : %s\n" "$OUT_DIR"
printf " finished   : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================================="

[[ "$fail" -ne 0 ]] && exit 1 || exit 0
