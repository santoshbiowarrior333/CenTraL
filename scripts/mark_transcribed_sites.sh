#!/usr/bin/env bash
# Mark distinct transcription sites per barcode.
#
# Built on top of step 5 — for each barcode it takes the per-read 1bp BED
# (barcodeNN.readstart.bed.gz, one row per primary read at its leftmost coord)
# and collapses to unique (chr, start, strand) tuples. The output is a
# presence/absence "transcribed sites" map — every position where at least
# one molecule was captured, with PCR amplification depth removed.
#
# For amplicon long-read data this is the right way to ask "where is
# transcription happening?" without `samtools markdup` over-collapsing
# legitimately-distinct molecules that happen to share a primer-defined
# 5' coordinate.
#
# Two files per barcode go into transcribed_sites/:
#   barcodeNN.sites.bed         BED6 — for intersections / karyoploteR ticks
#   barcodeNN.sites.bedgraph    chr/start/end/1 — drop-in for the existing
#                               karyoplot script (every bar has height 1, so
#                               the plot shows "where" not "how much")
#
# Plus a summary TSV with read-to-unique-site ratios per barcode.
#
# Activate your env first (samtools / gzip / awk / sort — all standard), then:
#     ./mark_transcribed_sites.sh -b /path/to/readstart_beds -t 8
#  or sbatch mark_transcribed_sites.sh   (edit BED_DIR below first)
#
# Positional shorthand:  ./mark_transcribed_sites.sh /path/to/readstart_beds

#SBATCH --job-name=mark_sites
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --time=0
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

# Default — edit for sbatch use, or override with -b
BED_DIR="/path/to/readstart_beds"
OUT_DIR=""
THREADS=""
FORCE=0
INCLUDE_UNCLASSIFIED=0

while getopts "b:o:t:fuh" opt; do
    case "$opt" in
        b) BED_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        f) FORCE=1 ;;
        u) INCLUDE_UNCLASSIFIED=1 ;;
        h) sed -n '2,30p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) exit 1 ;;
    esac
done

shift $((OPTIND - 1))
[[ $# -gt 0 && -d "$1" ]] && BED_DIR="$1"

[[ -z "$THREADS" ]] && THREADS="${SLURM_CPUS_PER_TASK:-4}"
BED_DIR="$(cd "$BED_DIR" && pwd)"
[[ -z "$OUT_DIR" ]] && OUT_DIR="$(dirname "$BED_DIR")/transcribed_sites"
mkdir -p "$OUT_DIR"

for tool in awk sort gzip; do
    if ! command -v "$tool" >/dev/null 2>&1; then
        echo "$tool not on PATH — that's surprising, check your env." >&2
        exit 127
    fi
done

echo "=================================================================="
echo " mark transcribed sites (unique per-position, strand-aware)"
echo "=================================================================="
echo " host         : $(hostname)"
echo " slurm job    : ${SLURM_JOB_ID:-(not in slurm)}"
echo " started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo " input  dir   : $BED_DIR"
echo " output dir   : $OUT_DIR"
echo " threads      : $THREADS   (used by sort)"
echo " force        : $FORCE"
echo " unclassified : $([[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && echo "included" || echo "skipped")"
echo "------------------------------------------------------------------"
echo

# Pick up the per-read BEDs from step 5. Same multi-glob style as the
# other scripts — barcode1 / barcode01 / barcode001 all match.
shopt -s nullglob
candidates=(
    "$BED_DIR"/barcode[0-9].readstart.bed.gz
    "$BED_DIR"/barcode[0-9][0-9].readstart.bed.gz
    "$BED_DIR"/barcode[0-9][0-9][0-9].readstart.bed.gz
)
[[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && candidates+=("$BED_DIR"/unclassified.readstart.bed.gz)
shopt -u nullglob

beds=()
for b in "${candidates[@]}"; do [[ -f "$b" ]] && beds+=("$b"); done

if [[ ${#beds[@]} -eq 0 ]]; then
    echo "No per-read BEDs in $BED_DIR (expected barcodeNN.readstart.bed.gz from step 5)." >&2
    exit 1
fi

echo "Found ${#beds[@]} BED(s) to process."
echo

summary="$OUT_DIR/transcribed_sites_summary.tsv"
printf "barcode\ttotal_reads\tunique_sites\tcollapse_ratio\n" > "$summary"

job_start=$SECONDS
done_count=0; skip_count=0; fail=0

for in_bed in "${beds[@]}"; do
    name="$(basename "$in_bed" .readstart.bed.gz)"
    sites_bed="$OUT_DIR/${name}.sites.bed"
    sites_bg="$OUT_DIR/${name}.sites.bedgraph"

    if [[ -s "$sites_bed" && -s "$sites_bg" && "$FORCE" -ne 1 ]]; then
        n_sites=$(wc -l < "$sites_bed")
        echo ">> $name: already done (${n_sites} sites) — skipping (use -f to overwrite)"
        skip_count=$((skip_count + 1)); echo; continue
    fi

    step_start=$SECONDS
    echo ">> $name"

    # Read count from the input BED — that's "reads with a primary alignment".
    in_reads=$(zcat -- "$in_bed" | wc -l)
    if [[ "$in_reads" -eq 0 ]]; then
        echo "   empty input — nothing to do for $name"
        skip_count=$((skip_count + 1)); echo; continue
    fi

    # Unique (chr, start, strand) — that's the "distinct transcription site"
    # definition. Same start on opposite strands = two sites (sense + antisense
    # of the same locus, biologically distinct).
    # awk fills BED6 with placeholder name/score, end = start + 1.
    if ! zcat -- "$in_bed" \
            | awk -v OFS='\t' '{print $1, $2, $3, ".", ".", $6}' \
            | sort --parallel="$THREADS" -k1,1 -k2,2n -k6,6 -u \
            > "$sites_bed"; then
        echo "   !! sites BED step FAILED for $name" >&2
        fail=$((fail + 1)); echo; continue
    fi

    # Drop-in bedgraph for karyoplot: every site gets height 1.
    # Sort key matches the BED (chr, start) so the bedgraph is coordinate-ordered.
    if ! awk -v OFS='\t' '{print $1, $2, $3, 1}' "$sites_bed" \
            | sort --parallel="$THREADS" -k1,1 -k2,2n -u \
            > "$sites_bg"; then
        echo "   !! sites bedgraph step FAILED for $name" >&2
        fail=$((fail + 1)); echo; continue
    fi

    n_sites=$(wc -l < "$sites_bed")
    ratio=$(awk -v r="$in_reads" -v s="$n_sites" 'BEGIN{ if(s>0) printf "%.2f", r/s; else print "NA" }')
    elapsed=$((SECONDS - step_start))
    printf "   reads %s   sites %s   reads-per-site %s   (%ds)\n\n" \
        "$in_reads" "$n_sites" "$ratio" "$elapsed"
    printf "%s\t%s\t%s\t%s\n" "$name" "$in_reads" "$n_sites" "$ratio" >> "$summary"
    done_count=$((done_count + 1))
done

total=$((SECONDS - job_start))
echo "=================================================================="
echo " Summary"
echo "------------------------------------------------------------------"
printf " processed  : %d\n" "$done_count"
printf " skipped    : %d\n" "$skip_count"
printf " failed     : %d\n" "$fail"
printf " total time : %dm %ds\n" $((total/60)) $((total%60))
printf " output dir : %s\n" "$OUT_DIR"
printf " summary    : %s\n" "$summary"
echo "------------------------------------------------------------------"
if [[ -s "$summary" ]]; then
    echo " per-barcode reads / unique sites / reads-per-site:"
    column -t -s $'\t' "$summary" 2>/dev/null | sed 's/^/   /' || sed 's/^/   /' "$summary"
fi
echo "=================================================================="

[[ "$fail" -ne 0 ]] && exit 1
exit 0
