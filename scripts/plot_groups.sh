#!/usr/bin/env bash
# Batch wrapper around karyoplot_bedgraph.R
#
# Loops over chromosomes and three hardcoded barcode groups
# (01-06, 07-12, 13-18), producing one karyoplot per (chromosome, group).
# Optional all-genome plot per group as well (-A).
#
# Assumes the standard layout:
#   <repo>/data_for_plots/       reference BEDs + chrom.sizes
#   <repo>/dcs_analysis_<ts>/    output of the main pipeline (passed via -d)
#   <repo>/plots/                where this script writes (auto-created)
#
# Examples:
#   ./scripts/plot_groups.sh -d dcs_analysis_20260524_130030
#   ./scripts/plot_groups.sh -d dcs_analysis_20260524_130030 -A
#   ./scripts/plot_groups.sh -d dcs_analysis_20260524_130030 -G     # only all-genome
#   ./scripts/plot_groups.sh -d dcs_analysis_20260524_130030 -C "chr1 chr2"
#
# To change the groups (e.g. you have 24 barcodes, or want non-contiguous
# treated/control sets), edit the GROUPS array below. Format per line:
#     "label:bc bc bc ..."   (space between barcodes, ":" between label and members)

#SBATCH --job-name=plot_groups
#SBATCH --cpus-per-task=4
#SBATCH --mem=0
#SBATCH --time=0
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
REPO_ROOT="$(cd "$SCRIPT_DIR/.." && pwd)"

# Defaults — edit if your layout differs
DATA_DIR="$REPO_ROOT/data_for_plots"
OUT_DIR="$REPO_ROOT/plots"
CHROM_SIZES="$DATA_DIR/rpe1_chrom.sizes"
REGIONS="$DATA_DIR/centromere_horAll.bed"
PRIMER_FWD="$DATA_DIR/forward_starts.bed"
PRIMER_REV="$DATA_DIR/reverse_starts.bed"
ASO_PLUS="$DATA_DIR/aso2_le3mm_plus.bed"
ASO_MINUS="$DATA_DIR/aso2_le3mm_minus.bed"

# Barcode groups. Edit these for your experiment.
GROUPS=(
    "b01-06:01 02 03 04 05 06"
    "b07-12:07 08 09 10 11 12"
    "b13-18:13 14 15 16 17 18"
)

DCS_DIR=""
CHROMS=()
DO_PER_CHR=1
DO_ALL=0

usage() {
    cat <<EOF
Usage: $(basename "$0") -d DCS_DIR [options]

  -d DIR    dcs_analysis_<timestamp> dir from the main pipeline (required)
  -o DIR    output dir (default: $OUT_DIR)
  -C "..."  space-separated chrom list (default: all in chrom.sizes)
  -A        ALSO produce a full-genome plot per group
  -G        ONLY produce the full-genome plot, skip per-chr
  -h        this help

To change groupings, edit the GROUPS array near the top of this script.
EOF
}

while getopts "d:o:C:AGh" opt; do
    case "$opt" in
        d) DCS_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        C) read -r -a CHROMS <<<"$OPTARG" ;;
        A) DO_ALL=1 ;;
        G) DO_ALL=1; DO_PER_CHR=0 ;;
        h) usage; exit 0 ;;
        *) usage >&2; exit 1 ;;
    esac
done

if [[ -z "$DCS_DIR" ]]; then
    echo "ERROR: -d DCS_DIR is required" >&2
    usage >&2
    exit 1
fi

DCS_DIR="$(cd "$DCS_DIR" && pwd)"
DCS_TSV="$DCS_DIR/dcs_counts.tsv"
BG_DIR="$DCS_DIR/readstart_beds"
KARYOPLOT="$SCRIPT_DIR/karyoplot_bedgraph.R"

for f in "$CHROM_SIZES" "$REGIONS" "$PRIMER_FWD" "$PRIMER_REV" \
         "$ASO_PLUS" "$ASO_MINUS" "$DCS_TSV" "$KARYOPLOT"; do
    [[ -f "$f" ]] || { echo "missing: $f" >&2; exit 1; }
done
[[ -d "$BG_DIR" ]] || { echo "missing dir: $BG_DIR" >&2; exit 1; }

if ! command -v Rscript >/dev/null 2>&1; then
    echo "Rscript not on PATH — activate your R env first." >&2
    exit 127
fi

if [[ ${#CHROMS[@]} -eq 0 ]]; then
    mapfile -t CHROMS < <(awk '{print $1}' "$CHROM_SIZES")
fi

mkdir -p "$OUT_DIR"

echo "=================================================================="
echo " batch karyoplot — per chr + barcode group"
echo "=================================================================="
echo " repo root  : $REPO_ROOT"
echo " dcs dir    : $DCS_DIR"
echo " out dir    : $OUT_DIR"
echo " chroms     : ${#CHROMS[@]}  (${CHROMS[0]} .. ${CHROMS[-1]})"
echo " groups     : ${#GROUPS[@]}"
for g in "${GROUPS[@]}"; do
    printf "              %-12s -> %s\n" "${g%%:*}" "${g#*:}"
done
echo " per-chr    : $([[ $DO_PER_CHR -eq 1 ]] && echo yes || echo no)"
echo " all-genome : $([[ $DO_ALL -eq 1 ]] && echo yes || echo no)"
echo "------------------------------------------------------------------"
echo

job_start=$SECONDS
total=0; ok=0; fail=0; skip=0

run_one() {
    local out="$1"; shift
    local chr="$1"; shift
    local zoom="$1"; shift
    local bgs=("$@")
    total=$((total + 1))

    if [[ ${#bgs[@]} -eq 0 ]]; then
        echo "  [$chr]  no input files — skipping"
        skip=$((skip + 1))
        return
    fi

    echo "  [$chr]  ${#bgs[@]} barcode(s)  ->  $out"
    if Rscript "$KARYOPLOT" \
        "$CHROM_SIZES" "$REGIONS" NA \
        "$out" "$chr" "$zoom" "$DCS_TSV" \
        "${bgs[@]}" \
        --primer-fwd "$PRIMER_FWD" --primer-rev "$PRIMER_REV" \
        --aso-plus   "$ASO_PLUS"   --aso-minus  "$ASO_MINUS" \
        > "${out}.log" 2>&1
    then
        ok=$((ok + 1))
    else
        echo "    !! FAILED — see ${out}.log" >&2
        fail=$((fail + 1))
    fi
}

if [[ $DO_PER_CHR -eq 1 ]]; then
    for g in "${GROUPS[@]}"; do
        label="${g%%:*}"
        members="${g#*:}"
        echo "## group $label"
        bgs=()
        for bc in $members; do
            f="$BG_DIR/barcode${bc}.startcount.bedgraph"
            [[ -f "$f" ]] && bgs+=("$f")
        done
        for chr in "${CHROMS[@]}"; do
            run_one "$OUT_DIR/${chr}_${label}" "$chr" "auto" "${bgs[@]}"
        done
        echo
    done
fi

if [[ $DO_ALL -eq 1 ]]; then
    echo "## all-genome plots"
    for g in "${GROUPS[@]}"; do
        label="${g%%:*}"
        members="${g#*:}"
        bgs=()
        for bc in $members; do
            f="$BG_DIR/barcode${bc}.startcount.bedgraph"
            [[ -f "$f" ]] && bgs+=("$f")
        done
        run_one "$OUT_DIR/all_${label}" "all" "NA" "${bgs[@]}"
    done
fi

elapsed=$((SECONDS - job_start))
echo
echo "=================================================================="
echo " Summary"
echo "------------------------------------------------------------------"
printf " attempted  : %d\n" "$total"
printf " ok         : %d\n" "$ok"
printf " skipped    : %d\n" "$skip"
printf " failed     : %d\n" "$fail"
printf " total time : %dm %ds\n" $((elapsed/60)) $((elapsed%60))
printf " out dir    : %s\n" "$OUT_DIR"
echo "=================================================================="

[[ "$fail" -ne 0 ]] && exit 1 || exit 0
