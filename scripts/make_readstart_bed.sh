#!/usr/bin/env bash
# Per-barcode read-start outputs from merged Nanopore BAMs. Two files per barcode:
#
#   1) <barcode>.readstart.bed.gz      BED6, ONE ENTRY PER PRIMARY READ
#         chrom  start  end  read_name  mapq  strand
#      Each interval is 1 bp wide at the read's leftmost genomic coordinate.
#      Good for inspecting individual reads or for downstream tools that want
#      per-read records.
#
#   2) <barcode>.startcount.bedgraph   COLLAPSED, ONE ENTRY PER UNIQUE POSITION
#         chrom  start  end  read_count
#      Standard bedgraph. Drop-in for karyoploteR (kpHeatmap / kpBars /
#      kpLines / kpArea), pyGenomeTracks, or anything that consumes bedgraph.
#
# Reads filtered with `samtools view -F 2308` — primary mapped only, no
# unmapped / secondary / supplementary. So each read contributes exactly one
# tick at the position where its alignment starts on the reference.
#
# Re-running is fast: if the per-read BED already exists, we skip the BAM
# extraction and just (re)derive the bedgraph from it.
#
# Activate your env first so samtools is on PATH, then:
#     ./make_readstart_bed.sh -b /path/to/merged_bam -t 8
#  or sbatch make_readstart_bed.sh   (edit BAM_DIR below first)

#SBATCH --job-name=read_starts
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
GZIP=1
INCLUDE_UNCLASSIFIED=0
FORCE=0

while getopts "b:o:t:Gufh" opt; do
    case "$opt" in
        b) BAM_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        G) GZIP=0 ;;          # write plain .bed instead of .bed.gz
        u) INCLUDE_UNCLASSIFIED=1 ;;
        f) FORCE=1 ;;
        h) sed -n '2,24p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) exit 1 ;;
    esac
done

shift $((OPTIND - 1))
[[ $# -gt 0 && -d "$1" ]] && BAM_DIR="$1"

[[ -z "$THREADS" ]] && THREADS="${SLURM_CPUS_PER_TASK:-4}"
BAM_DIR="$(cd "$BAM_DIR" && pwd)"
[[ -z "$OUT_DIR" ]] && OUT_DIR="$(dirname "$BAM_DIR")/readstart_beds"
mkdir -p "$OUT_DIR"

for t in samtools awk sort uniq; do
    command -v "$t" >/dev/null 2>&1 || { echo "$t not on PATH — activate your env." >&2; exit 127; }
done

echo "=================================================================="
echo " read-start outputs (per-read BED + collapsed bedgraph)"
echo "=================================================================="
echo " host         : $(hostname)"
echo " slurm job    : ${SLURM_JOB_ID:-(not in slurm)}"
echo " started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo " samtools     : $(samtools --version | head -n1)"
echo " bam dir      : $BAM_DIR"
echo " output dir   : $OUT_DIR"
echo " threads      : $THREADS"
echo " filter       : -F 2308  (skip unmapped + secondary + supplementary)"
echo " gzip BED     : $GZIP"
echo " unclassified : $([[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && echo "included" || echo "skipped")"
echo "------------------------------------------------------------------"
echo

shopt -s nullglob
candidates=("$BAM_DIR"/barcode[0-9][0-9].bam)
[[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && candidates+=("$BAM_DIR"/unclassified.bam)
shopt -u nullglob
bams=()
for b in "${candidates[@]}"; do [[ -f "$b" ]] && bams+=("$b"); done

if [[ ${#bams[@]} -eq 0 ]]; then
    echo "No barcode BAMs in $BAM_DIR (expected barcodeNN.bam / unclassified.bam)." >&2
    exit 1
fi

echo "Found ${#bams[@]} BAM(s) to process."
echo

job_start=$SECONDS
done_count=0; skip_count=0; fail=0

# Helper: read a BED file (gzipped or not) to stdout
cat_bed() {
    if [[ "$1" == *.gz ]]; then zcat "$1"; else cat "$1"; fi
}

for b in "${bams[@]}"; do
    name="$(basename "$b" .bam)"
    out_bed="$OUT_DIR/${name}.readstart.bed"
    [[ "$GZIP" -eq 1 ]] && out_bed="${out_bed}.gz"
    out_bg="$OUT_DIR/${name}.startcount.bedgraph"

    if [[ -s "$out_bed" && -s "$out_bg" && "$FORCE" -ne 1 ]]; then
        echo ">> $name: both outputs already exist — skipping (use -f to overwrite)"
        skip_count=$((skip_count + 1)); echo; continue
    fi

    step_start=$SECONDS
    echo ">> $name"

    # --- 1. per-read BED ---
    if [[ -s "$out_bed" && "$FORCE" -ne 1 ]]; then
        echo "   per-read BED already exists — keeping it"
        n=$(cat_bed "$out_bed" | wc -l)
    else
        echo -n "   extracting read starts from BAM..."
        # samtools -F 2308 = primary mapped only. BAM POS is 1-based; BED start is 0-based.
        # bit 16 of FLAG => reverse strand.
        if [[ "$GZIP" -eq 1 ]]; then
            samtools view -@ "$THREADS" -F 2308 "$b" \
              | awk 'BEGIN{OFS="\t"}
                     { strand = (int($2/16) % 2 == 1) ? "-" : "+"
                       print $3, $4-1, $4, $1, $5, strand }' \
              | gzip -c > "$out_bed"
        else
            samtools view -@ "$THREADS" -F 2308 "$b" \
              | awk 'BEGIN{OFS="\t"}
                     { strand = (int($2/16) % 2 == 1) ? "-" : "+"
                       print $3, $4-1, $4, $1, $5, strand }' \
              > "$out_bed"
        fi
        n=$(cat_bed "$out_bed" | wc -l)
        size=$(du -h "$out_bed" | cut -f1)
        echo " done — $n entries, $size"
    fi

    # --- 2. collapsed bedgraph ---
    if [[ -s "$out_bg" && "$FORCE" -ne 1 ]]; then
        echo "   bedgraph already exists — keeping it"
        npos=$(wc -l < "$out_bg")
    else
        echo -n "   collapsing to bedgraph..."
        cat_bed "$out_bed" \
          | awk 'BEGIN{OFS="\t"} {print $1, $2, $3}' \
          | sort -k1,1 -k2,2n \
          | uniq -c \
          | awk 'BEGIN{OFS="\t"} {print $2, $3, $4, $1}' \
          > "$out_bg"
        npos=$(wc -l < "$out_bg")
        size=$(du -h "$out_bg" | cut -f1)
        echo " done — $npos unique positions, $size"
    fi

    elapsed=$((SECONDS - step_start))
    printf "   %dm %ds total\n\n" $((elapsed/60)) $((elapsed%60))
    done_count=$((done_count + 1))
done

total=$((SECONDS - job_start))
echo "=================================================================="
echo " Summary"
echo "------------------------------------------------------------------"
printf " processed   : %d\n" "$done_count"
printf " skipped     : %d\n" "$skip_count"
printf " failed      : %d\n" "$fail"
printf " total time  : %dm %ds\n" $((total/60)) $((total%60))
printf " output dir  : %s\n" "$OUT_DIR"
printf " files       : <barcode>.readstart.bed[.gz]  (per-read BED6)\n"
printf "               <barcode>.startcount.bedgraph (chr/start/end/count, for karyoploteR)\n"
printf " finished    : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================================="
