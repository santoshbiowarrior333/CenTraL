#!/usr/bin/env bash
# Generate per-barcode bigwig coverage tracks from merged Nanopore BAMs,
# scaled by the DCS spike-in factor from count_dcs_spikein.sh.
#
# Reads dcs_counts.tsv (column 'scale_factor') and runs deepTools bamCoverage
# on every barcode that has a usable scale factor. Skips barcodes marked NA.
#
# Activate your env first so deepTools (bamCoverage) and samtools are on PATH:
#     conda activate python_collection
# Then run interactively or sbatch.
#
# Interactive: ./make_normalized_bigwigs.sh -b /path/to/merged_bam -t 8
# SLURM:       sbatch make_normalized_bigwigs.sh  (edit BAM_DIR below first)
#
# Common usage patterns for a diploid (hap-resolved) reference:
#   # all reads, multimappers included — for whole-genome visualization
#   ./make_normalized_bigwigs.sh -b merged_bam -q 0  -o bw_all
#   # high-MAPQ only — for hap-resolved/allele-specific signal
#   ./make_normalized_bigwigs.sh -b merged_bam -q 20 -o bw_hq

#SBATCH --job-name=make_bw
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=0
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

# Default paths — edit BAM_DIR for sbatch use, or override with -b
BAM_DIR="/path/to/merged_bam"

# CLI defaults (overridable with flags)
OUT_DIR=""
TSV=""
THREADS=""
BIN=50
MIN_MAPQ=0
FORCE=0

while getopts "b:o:c:t:s:q:fh" opt; do
    case "$opt" in
        b) BAM_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        c) TSV="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        s) BIN="$OPTARG" ;;
        q) MIN_MAPQ="$OPTARG" ;;
        f) FORCE=1 ;;
        h) sed -n '2,20p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) exit 1 ;;
    esac
done

shift $((OPTIND - 1))
[[ $# -gt 0 && -d "$1" ]] && BAM_DIR="$1"

[[ -z "$THREADS" ]] && THREADS="${SLURM_CPUS_PER_TASK:-4}"
BAM_DIR="$(cd "$BAM_DIR" && pwd)"
[[ -z "$TSV"     ]] && TSV="$BAM_DIR/dcs_counts.tsv"
[[ -z "$OUT_DIR" ]] && OUT_DIR="$(dirname "$BAM_DIR")/bw_mapq${MIN_MAPQ}"
mkdir -p "$OUT_DIR"

# Required tools
for t in bamCoverage samtools awk; do
    if ! command -v "$t" >/dev/null 2>&1; then
        echo "$t not on PATH — activate the right conda env first." >&2
        exit 127
    fi
done

if [[ ! -f "$TSV" ]]; then
    echo "DCS counts TSV not found: $TSV" >&2
    echo "Run count_dcs_spikein.sh first to produce it, or pass an alternate path with -c." >&2
    exit 1
fi

# Banner
echo "=================================================================="
echo " bigwig generation (DCS spike-in normalized)"
echo "=================================================================="
echo " host          : $(hostname)"
echo " slurm job     : ${SLURM_JOB_ID:-(not in slurm)}"
echo " started       : $(date '+%Y-%m-%d %H:%M:%S')"
echo " bamCoverage   : $(bamCoverage --version 2>&1 | tail -n1)"
echo " bam dir       : $BAM_DIR"
echo " dcs tsv       : $TSV"
echo " output dir    : $OUT_DIR"
echo " bin size      : $BIN bp"
if [[ "$MIN_MAPQ" -ge 20 ]]; then
    mq_note="(hap-resolved / high-MAPQ only — multimappers excluded)"
elif [[ "$MIN_MAPQ" -gt 0 ]]; then
    mq_note="(MAPQ filter active)"
else
    mq_note="(all reads — multimappers included)"
fi
echo " min MAPQ      : $MIN_MAPQ  $mq_note"
echo " threads       : $THREADS"
echo " force         : $FORCE"
echo "------------------------------------------------------------------"
echo

# Slurp the DCS table — barcode name and scale factor are all we need.
declare -A SCALES
declare -A DCS_COUNTS
order=()
while IFS=$'\t' read -r bc _ _ _ dcs sf; do
    [[ "$bc" == "barcode" || -z "$bc" ]] && continue
    SCALES["$bc"]="$sf"
    DCS_COUNTS["$bc"]="$dcs"
    order+=("$bc")
done < "$TSV"

if [[ ${#order[@]} -eq 0 ]]; then
    echo "No barcodes parsed from $TSV." >&2
    exit 1
fi

# Pre-flight: print the plan so we can sanity-check before the long bit
echo "Plan (${#order[@]} barcode(s) in TSV):"
for bc in "${order[@]}"; do
    sf="${SCALES[$bc]}"
    dcs="${DCS_COUNTS[$bc]}"
    bam="$BAM_DIR/${bc}.bam"
    if [[ "$sf" == "NA" ]]; then
        printf "  %-15s  skip — no DCS reads (scale=NA)\n" "$bc"
    elif [[ ! -f "$bam" ]]; then
        printf "  %-15s  skip — BAM not found in %s\n" "$bc" "$BAM_DIR"
    else
        printf "  %-15s  scale=%-10s dcs=%s\n" "$bc" "$sf" "$dcs"
    fi
done
echo

# Run per-barcode
job_start=$SECONDS
done_count=0
skip_count=0
fail=0

for bc in "${order[@]}"; do
    sf="${SCALES[$bc]}"
    bam="$BAM_DIR/${bc}.bam"
    bw="$OUT_DIR/${bc}.bw"

    if [[ "$sf" == "NA" ]]; then
        echo ">> $bc: scale=NA, skipping"
        skip_count=$((skip_count + 1))
        echo
        continue
    fi
    if [[ ! -f "$bam" ]]; then
        echo ">> $bc: BAM not found ($bam), skipping"
        skip_count=$((skip_count + 1))
        echo
        continue
    fi
    if [[ -s "$bw" && "$FORCE" -ne 1 ]]; then
        echo ">> $bc: bigwig already exists ($(du -h "$bw" | cut -f1)), use -f to overwrite — skipping"
        skip_count=$((skip_count + 1))
        echo
        continue
    fi

    step_start=$SECONDS
    echo ">> $bc"
    echo "   scale factor : $sf"
    echo "   --> $bw"
    echo -n "   running bamCoverage..."

    cmd=(bamCoverage
         -b "$bam"
         -o "$bw"
         --scaleFactor "$sf"
         --binSize "$BIN"
         -p "$THREADS")
    [[ "$MIN_MAPQ" -gt 0 ]] && cmd+=(--minMappingQuality "$MIN_MAPQ")

    log="$OUT_DIR/${bc}.bamCoverage.log"
    if ! "${cmd[@]}" >"$log" 2>&1; then
        echo " FAILED (see $log)"
        fail=$((fail + 1))
        echo
        continue
    fi

    elapsed=$((SECONDS - step_start))
    size=$(du -h "$bw" | cut -f1)
    echo " done"
    printf "   size %s, %dm %ds\n\n" "$size" $((elapsed/60)) $((elapsed%60))
    done_count=$((done_count + 1))
done

total_elapsed=$((SECONDS - job_start))
echo "=================================================================="
echo " Summary"
echo "------------------------------------------------------------------"
printf " bigwigs generated : %d\n" "$done_count"
printf " skipped           : %d\n" "$skip_count"
printf " failed            : %d\n" "$fail"
printf " total time        : %dm %ds\n" $((total_elapsed/60)) $((total_elapsed%60))
printf " output dir        : %s\n" "$OUT_DIR"
printf " bamCoverage logs  : %s/*.bamCoverage.log\n" "$OUT_DIR"
printf " finished          : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================================="

[[ "$fail" -ne 0 ]] && exit 1 || exit 0
