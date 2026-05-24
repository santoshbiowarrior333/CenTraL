#!/usr/bin/env bash
# Full Nanopore PCR-amplicon pipeline — one command, starts from raw bam_pass/.
#
# Pass the bam_pass/ directory (the standard MinKNOW/Dorado output with
# barcodeNN/ subfolders of chunked BAMs) and this runs:
#
#     Step 1 — merge chunked per-barcode BAMs    (-> bam_pass/merged_bam/)
#     Step 2 — count DCS spike-in                (-> dcs_analysis_TS/dcs_counts.tsv)
#     Step 3 — primary-only filter (-F 2308)     (-> dcs_analysis_TS/primary_bams/)
#     Step 4 — normalized bigwigs                (-> dcs_analysis_TS/bw/)
#     Step 5 — 1bp BEDs + bedgraphs              (-> dcs_analysis_TS/readstart_beds/)
#     Step 6 — DCS-normalization QC plot         (-> dcs_analysis_TS/dcs_normalization_qc.{png,pdf})
#
# The karyoplot step (step 7) is intentionally NOT included here — it needs
# user-specific chrom.sizes + HOR/centromere BEDs that vary per project.
# Run it separately with scripts/karyoplot_bedgraph.R after this finishes.
# See the README for examples.
#
# Activate your env first so the tools are on PATH:
#     conda activate python_collection      # or your env that has the deps
#     module load dorado                    # if dorado lives in a module
# Then:
#     ./run_dcs_workflow.sh /path/to/bam_pass
#  or sbatch run_dcs_workflow.sh /path/to/bam_pass

#SBATCH --job-name=dcs_workflow
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=0
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

# Default — edit for sbatch use, or override with -i
BAM_PASS_DIR="/path/to/bam_pass"

OUT_DIR=""
THREADS=""
BIN_SIZE=50
FORCE=0
INCLUDE_UNCLASSIFIED=0

while getopts "i:o:t:s:fuh" opt; do
    case "$opt" in
        i) BAM_PASS_DIR="$OPTARG" ;;
        o) OUT_DIR="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        s) BIN_SIZE="$OPTARG" ;;
        f) FORCE=1 ;;
        u) INCLUDE_UNCLASSIFIED=1 ;;
        h) sed -n '2,28p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) exit 1 ;;
    esac
done

# Allow positional first arg as shorthand for -i
shift $((OPTIND - 1))
[[ $# -gt 0 && -d "$1" ]] && BAM_PASS_DIR="$1"

[[ -z "$THREADS" ]] && THREADS="${SLURM_CPUS_PER_TASK:-4}"

BAM_PASS_DIR="$(cd "$BAM_PASS_DIR" && pwd)"
SCRIPT_DIR="$(dirname "$(readlink -f "$0")")"
SUBSCRIPTS="$SCRIPT_DIR/scripts"
DCS_REF="$SUBSCRIPTS/DCS_Lambda_3.6kb.fa"
MERGED_BAM_DIR="$BAM_PASS_DIR/merged_bam"

if [[ -z "$OUT_DIR" ]]; then
    ts="$(date '+%Y%m%d_%H%M%S')"
    OUT_DIR="$(pwd)/dcs_analysis_${ts}"
fi
mkdir -p "$OUT_DIR"

# Preflight — fail fast before SLURM allocates anything
missing=()
for t in samtools bamCoverage awk python3 gzip; do
    command -v "$t" >/dev/null 2>&1 || missing+=("$t")
done
for pymod in matplotlib numpy; do
    python3 -c "import $pymod" 2>/dev/null || missing+=("python:$pymod")
done
if ! command -v dorado >/dev/null 2>&1 && ! command -v minimap2 >/dev/null 2>&1; then
    missing+=("dorado-or-minimap2")
fi
# Make sure all the sub-scripts are where we expect them
for s in barcode_merge_nanopore.sh count_dcs_spikein.sh \
         filter_primary_bams.sh make_normalized_bigwigs.sh \
         make_readstart_bed.sh plot_normalization_qc.py; do
    [[ -f "$SUBSCRIPTS/$s" ]] || missing+=("scripts/$s")
done

if [[ ${#missing[@]} -gt 0 ]]; then
    echo "Preflight FAILED — missing:" >&2
    for m in "${missing[@]}"; do echo "  - $m" >&2; done
    echo "Activate the right env / load the right modules and try again." >&2
    echo "On cluster1: conda activate python_collection" >&2
    exit 127
fi

# Tee everything from here on into run.log inside the output dir
LOG="$OUT_DIR/run.log"
exec > >(tee -a "$LOG") 2>&1

echo "=================================================================="
echo " Nanopore PCR-amplicon pipeline — full run"
echo "=================================================================="
echo " host         : $(hostname)"
echo " slurm job    : ${SLURM_JOB_ID:-(not in slurm)}"
echo " started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo " bam_pass     : $BAM_PASS_DIR"
echo " merged_bam   : $MERGED_BAM_DIR"
echo " output dir   : $OUT_DIR"
echo " threads      : $THREADS"
echo " bin size     : $BIN_SIZE"
echo " unclassified : $([[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && echo "included" || echo "skipped after step 1")"
echo "=================================================================="

job_start=$SECONDS

# ---- Step 1: merge per-barcode chunked BAMs -> bam_pass/merged_bam/ ----
echo
echo "+++ Step 1/6 — merge chunked per-barcode BAMs +++"
merge_flags=(-i "$BAM_PASS_DIR" -t "$THREADS")
[[ "$FORCE" -eq 1 ]] && merge_flags+=(-f)
bash "$SUBSCRIPTS/barcode_merge_nanopore.sh" "${merge_flags[@]}"

# ---- Step 2: count DCS spike-in -> OUT_DIR/dcs_counts.tsv ----
echo
echo "+++ Step 2/6 — count DCS spike-in +++"
TSV="$OUT_DIR/dcs_counts.tsv"
count_flags=(-b "$MERGED_BAM_DIR" -r "$DCS_REF" -o "$TSV" -t "$THREADS")
[[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && count_flags+=(-u)
bash "$SUBSCRIPTS/count_dcs_spikein.sh" "${count_flags[@]}"

# ---- Step 3: primary-only filter -> OUT_DIR/primary_bams/ ----
echo
echo "+++ Step 3/6 — filter to primary mapped reads (-F 2308) +++"
PRIMARY_DIR="$OUT_DIR/primary_bams"
fp_flags=(-b "$MERGED_BAM_DIR" -o "$PRIMARY_DIR" -t "$THREADS")
[[ "$FORCE" -eq 1 ]] && fp_flags+=(-f)
[[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && fp_flags+=(-u)
bash "$SUBSCRIPTS/filter_primary_bams.sh" "${fp_flags[@]}"

# ---- Step 4: normalized bigwigs from primary_bams ----
echo
echo "+++ Step 4/6 — normalized bigwigs +++"
BW_DIR="$OUT_DIR/bw"
bw_flags=(-b "$PRIMARY_DIR" -o "$BW_DIR" -c "$TSV" -t "$THREADS" -q 0 -s "$BIN_SIZE")
[[ "$FORCE" -eq 1 ]] && bw_flags+=(-f)
bash "$SUBSCRIPTS/make_normalized_bigwigs.sh" "${bw_flags[@]}"

# ---- Step 5: read-start BEDs + bedgraphs ----
echo
echo "+++ Step 5/6 — read-start BEDs + bedgraphs (1 bp per primary read) +++"
BED_DIR="$OUT_DIR/readstart_beds"
bed_flags=(-b "$PRIMARY_DIR" -o "$BED_DIR" -t "$THREADS")
[[ "$FORCE" -eq 1 ]] && bed_flags+=(-f)
[[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && bed_flags+=(-u)
bash "$SUBSCRIPTS/make_readstart_bed.sh" "${bed_flags[@]}"

# ---- Step 6: QC plot ----
echo
echo "+++ Step 6/6 — DCS normalization QC plot +++"
python3 "$SUBSCRIPTS/plot_normalization_qc.py" \
    -c "$TSV" \
    -o "$OUT_DIR/dcs_normalization_qc"

# ---- summary ----
total_elapsed=$((SECONDS - job_start))
echo
echo "=================================================================="
echo " ALL DONE"
echo "------------------------------------------------------------------"
printf " total time : %dh %dm %ds\n" \
    $((total_elapsed/3600)) $(((total_elapsed%3600)/60)) $((total_elapsed%60))
printf " finished   : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
printf " output dir : %s\n" "$OUT_DIR"
echo "------------------------------------------------------------------"
echo " key outputs:"
echo "   merged_bam/ (in source tree, $BAM_PASS_DIR/merged_bam)"
echo "   $OUT_DIR/dcs_counts.tsv"
echo "   $OUT_DIR/primary_bams/*.bam"
echo "   $OUT_DIR/bw/*.bw                              (normalized for IGV)"
echo "   $OUT_DIR/readstart_beds/*.startcount.bedgraph (for karyoploteR)"
echo "   $OUT_DIR/dcs_normalization_qc.{png,pdf}"
echo
echo " Next: per-chromosome karyoplots — see scripts/karyoplot_bedgraph.R"
echo "       (needs your project's chrom.sizes + HOR/centromere BEDs)"
echo "=================================================================="
