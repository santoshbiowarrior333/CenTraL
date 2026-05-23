#!/usr/bin/env bash
# Count DCS (ONT lambda spike-in) reads in merged Nanopore BAMs and emit a
# tidy TSV with a scaling factor you can feed straight into bamCoverage.
#
# For each barcodeNN.bam in the bam dir, this:
#   1) pulls reads that did NOT map to your target genome (the unmapped pool)
#   2) re-aligns them to the DCS reference (3.6 kb lambda fragment)
#   3) counts how many came back mapped — that's your DCS count
#   4) writes dcs_counts.tsv: barcode | total | target_mapped | unmapped | dcs | scale
#
# Scaling factor:  scale = min(dcs_across_samples) / dcs_this_sample
# i.e. the sample with the fewest DCS reads gets scale=1.0, others get scaled DOWN.
# Use it directly with deepTools: bamCoverage --scaleFactor <scale> ...
#
# The DCS reference is embedded in this script — if the FASTA isn't found at
# the expected path, it's written from the heredoc at the bottom and indexed
# automatically with samtools faidx. Nothing to download or stage.
#
# Activate your env first so samtools + minimap2 are on PATH:
#     conda activate python_collection
# Then run interactively or submit with sbatch.
#
# Interactive: ./count_dcs_spikein.sh -b /path/to/merged_bam -t 8
# SLURM:       sbatch count_dcs_spikein.sh    (edit BAM_DIR below first)

#SBATCH --job-name=dcs_count
#SBATCH --cpus-per-task=8
#SBATCH --mem=0
#SBATCH --time=0
#SBATCH --output=slurm-%j.out
#SBATCH --error=slurm-%j.err

set -euo pipefail

# Default paths — edit BAM_DIR for sbatch use, or override at runtime with -b
BAM_DIR="/path/to/merged_bam"
DCS_REF="$(dirname "$(readlink -f "$0")")/DCS_Lambda_3.6kb.fa"

OUT_TSV=""
THREADS=""
INCLUDE_UNCLASSIFIED=0

while getopts "b:r:o:t:uh" opt; do
    case "$opt" in
        b) BAM_DIR="$OPTARG" ;;
        r) DCS_REF="$OPTARG" ;;
        o) OUT_TSV="$OPTARG" ;;
        t) THREADS="$OPTARG" ;;
        u) INCLUDE_UNCLASSIFIED=1 ;;
        h) sed -n '2,25p' "$0" | sed 's/^# \{0,1\}//'; exit 0 ;;
        *) exit 1 ;;
    esac
done

shift $((OPTIND - 1))
# Allow positional arg as shorthand for -b
[[ $# -gt 0 && -d "$1" ]] && BAM_DIR="$1"

[[ -z "$THREADS" ]] && THREADS="${SLURM_CPUS_PER_TASK:-4}"
BAM_DIR="$(cd "$BAM_DIR" && pwd)"
[[ -z "$OUT_TSV" ]] && OUT_TSV="$BAM_DIR/dcs_counts.tsv"

# Required tools — samtools + awk are mandatory; aligner is dorado OR minimap2.
for t in samtools awk; do
    if ! command -v "$t" >/dev/null 2>&1; then
        echo "$t not on PATH — activate the right conda env first." >&2
        exit 127
    fi
done

if command -v dorado >/dev/null 2>&1; then
    ALIGNER="dorado"
    ALIGNER_VER="$(dorado --version 2>&1 | head -n1)"
elif command -v minimap2 >/dev/null 2>&1; then
    ALIGNER="minimap2"
    ALIGNER_VER="$(minimap2 --version)"
else
    echo "Neither dorado nor minimap2 on PATH — load one of them and try again." >&2
    exit 127
fi

# ---- Bootstrap the DCS reference if it's not on disk yet ----
# Heredoc with the sequence is at the bottom of this file; we write from a
# delimited block so the script stays self-contained.
write_embedded_dcs_fa() {
    local out="$1"
    # Find the line range between the BEGIN/END markers and emit the FASTA
    awk '/^# >>> EMBEDDED_DCS_FASTA_BEGIN >>>$/{flag=1;next}
         /^# <<< EMBEDDED_DCS_FASTA_END <<<$/{flag=0}
         flag{sub(/^# ?/,""); print}' "$0" > "$out"
}

if [[ ! -f "$DCS_REF" ]]; then
    echo "DCS reference not found at $DCS_REF — writing embedded copy"
    if ! write_embedded_dcs_fa "$DCS_REF" 2>/dev/null; then
        # Target dir wasn't writable — fall back to a temp file
        DCS_REF="$(mktemp --suffix=.DCS_Lambda_3.6kb.fa)"
        echo "  (target not writable, using $DCS_REF instead)"
        write_embedded_dcs_fa "$DCS_REF"
    fi
fi

# Make sure samtools faidx index exists
if [[ ! -f "${DCS_REF}.fai" ]]; then
    samtools faidx "$DCS_REF"
fi

# Per-barcode minimap2 logs go here so the main stdout stays readable.
LOG_DIR="$(dirname "$OUT_TSV")/dcs_logs"
mkdir -p "$LOG_DIR"

echo "=================================================================="
echo " DCS spike-in counter"
echo "=================================================================="
echo " host         : $(hostname)"
echo " slurm job    : ${SLURM_JOB_ID:-(not in slurm)}"
echo " started      : $(date '+%Y-%m-%d %H:%M:%S')"
echo " samtools     : $(samtools --version | head -n1)"
echo " aligner      : $ALIGNER  ($ALIGNER_VER)"
echo " bam dir      : $BAM_DIR"
echo " dcs ref      : $DCS_REF ($(awk '/^>/{next}{t+=length($0)}END{print t}' "$DCS_REF") bp)"
echo " output tsv   : $OUT_TSV"
echo " logs dir     : $LOG_DIR"
echo " threads      : $THREADS"
echo "------------------------------------------------------------------"
echo

# Find barcode BAMs (and unclassified) in the directory.
shopt -s nullglob
candidates=("$BAM_DIR"/barcode[0-9][0-9].bam)
[[ "$INCLUDE_UNCLASSIFIED" -eq 1 ]] && candidates+=("$BAM_DIR"/unclassified.bam)
shopt -u nullglob

bams=()
for b in "${candidates[@]}"; do
    [[ -f "$b" ]] && bams+=("$b")
done

if [[ ${#bams[@]} -eq 0 ]]; then
    echo "No barcode BAMs found in $BAM_DIR (expected barcodeNN.bam / unclassified.bam)." >&2
    exit 1
fi

echo "Found ${#bams[@]} BAM(s) to scan."
echo

RAW=$(mktemp)
trap 'rm -f "$RAW"' EXIT
printf 'barcode\ttotal_reads\ttarget_mapped\tunmapped\tdcs_mapped\n' > "$RAW"

job_start=$SECONDS

for b in "${bams[@]}"; do
    name="$(basename "$b" .bam)"
    step_start=$SECONDS
    echo ">> $name"

    total=$(samtools view -@ "$THREADS" -c "$b")
    mapped=$(samtools view -@ "$THREADS" -c -F 4 "$b")
    unmapped=$((total - mapped))

    printf "   total reads     : %s\n" "$total"
    printf "   mapped (target) : %s\n" "$mapped"
    printf "   unmapped pool   : %s\n" "$unmapped"

    if [[ "$unmapped" -eq 0 ]]; then
        echo "   (no unmapped reads — nothing to align to DCS)"
        dcs=0
    else
        echo -n "   aligning unmapped reads to DCS (using $ALIGNER)..."
        log="$LOG_DIR/${name}.${ALIGNER}.log"

        # The merged BAM is already mapped — we extract just the unmapped
        # fraction and feed it to the aligner. dorado is happy with a BAM
        # input (its native format and it preserves any aux tags); minimap2
        # needs fastq. Stage either to a temp file and rm at the end.
        if [[ "$ALIGNER" == "dorado" ]]; then
            tmp="$(mktemp --suffix=.unmapped.bam)"
            samtools view -@ "$THREADS" -b -f 4 "$b" -o "$tmp" 2>>"$log"
            if [[ ! -s "$tmp" ]]; then
                dcs=0
            else
                # `-F 2308` = primary mapped only (no secondary / supplementary)
                dcs=$(dorado aligner --threads "$THREADS" "$DCS_REF" "$tmp" 2>>"$log" \
                      | samtools view -@ "$THREADS" -c -F 2308 -)
            fi
        else
            tmp="$(mktemp --suffix=.unmapped.fq)"
            samtools fastq -@ "$THREADS" -f 4 "$b" 2>>"$log" > "$tmp"
            if [[ ! -s "$tmp" ]]; then
                dcs=0
            else
                dcs=$(minimap2 -ax map-ont -t "$THREADS" "$DCS_REF" "$tmp" 2>>"$log" \
                      | samtools view -@ "$THREADS" -c -F 2308 -)
            fi
        fi
        rm -f "$tmp"
        echo " done"
    fi

    elapsed=$((SECONDS - step_start))
    if [[ "$total" -gt 0 ]]; then
        pct=$(awk -v d="$dcs" -v t="$total" 'BEGIN{printf "%.3f", 100*d/t}')
    else
        pct="0"
    fi
    printf "   DCS reads       : %s  (%s%% of total)\n" "$dcs" "$pct"
    printf "   step time       : %dm %ds\n\n" $((elapsed/60)) $((elapsed%60))

    printf '%s\t%s\t%s\t%s\t%s\n' "$name" "$total" "$mapped" "$unmapped" "$dcs" >> "$RAW"
done

echo "Writing final table with scaling factors..."
awk -F'\t' '
NR==1 { header=$0; next }
{
    rows[NR]=$0
    dcs[NR]=$5
    if (dcs[NR] > 0 && (min == 0 || dcs[NR] < min)) min = dcs[NR]
}
END {
    if (min == 0) min = 1
    print header "\tscale_factor"
    for (i = 2; i <= NR; i++) {
        if (dcs[i] == 0) sf = "NA"
        else               sf = sprintf("%.6f", min / dcs[i])
        print rows[i] "\t" sf
    }
}' "$RAW" > "$OUT_TSV"

echo
echo "----- $OUT_TSV -----"
column -t -s $'\t' "$OUT_TSV"
echo

total_elapsed=$((SECONDS - job_start))
echo "=================================================================="
echo " Summary"
echo "------------------------------------------------------------------"
printf " barcodes scanned : %d\n" "${#bams[@]}"
printf " total time       : %dm %ds\n" $((total_elapsed/60)) $((total_elapsed%60))
printf " output table     : %s\n" "$OUT_TSV"
printf " aligner logs     : %s/\n" "$LOG_DIR"
printf " finished         : %s\n" "$(date '+%Y-%m-%d %H:%M:%S')"
echo "=================================================================="
echo
echo "Next step — use the scale_factor column with deepTools bamCoverage:"
echo "  bamCoverage -b merged_bam/barcode01.bam -o bw/barcode01.bw \\"
echo "      --scaleFactor <scale_from_tsv> --binSize 50 -p $THREADS"

exit 0

# Everything below this point is data, not code. The embedded DCS FASTA
# lives between the BEGIN/END markers and is extracted by write_embedded_dcs_fa().
# Lines start with '# ' so the rest of the file is still pure bash.
#
# >>> EMBEDDED_DCS_FASTA_BEGIN >>>
# >DCS_Lambda_3.6kb
# GCCATCAGATTGTGTTTGTTAGTCGCTGCCATCAGATTGTGTTTGTTAGTCGCTTTTTTTTTTTGGAATTTTTTTTTTGGAATTTTTTTTTTGCGCTAACAACCTCCTGCCGTTTTGCCCGTGCATATCGGTCACGAACAAATCTGATTACTAAACACAGTAGCCTGGATTTGTTCTATCAGTAATCGACCTTATTCCTAATTAAATAGAGCAAATCCCCTTATTGGGGGTAAGACATGAAGATGCCAGAAAAACATGACCTGTTGGCCGCCATTCTCGCGGCAAAGGAACAAGGCATCGGGGCAATCCTTGCGTTTGCAATGGCGTACCTTCGCGGCAGATATAATGGCGGTGCGTTTACAAAAACAGTAATCGACGCAACGATGTGCGCCATTATCGCCTAGTTCATTCGTGACCTTCTCGACTTCGCCGGACTAAGTAGCAATCTCGCTTATATAACGAGCGTGTTTATCGGCTACATCGGTACTGACTCGATTGGTTCGCTTATCAAACGCTTCGCTGCTAAAAAAGCCGGAGTAGAAGATGGTAGAAATCAATAATCAACGTAAGGCGTTCCTCGATATGCTGGCGTGGTCGGAGGGAACTGATAACGGACGTCAGAAAACCAGAAATCATGGTTATGACGTCATTGTAGGCGGAGAGCTATTTACTGATTACTCCGATCACCCTCGCAAACTTGTCACGCTAAACCCAAAACTCAAATCAACAGGCGCCGGACGCTACCAGCTTCTTTCCCGTTGGTGGGATGCCTACCGCAAGCAGCTTGGCCTGAAAGACTTCTCTCCGAAAAGTCAGGACGCTGTGGCATTGCAGCAGATTAAGGAGCGTGGCGCTTTACCTATGATTGATCGTGGTGATATCCGTCAGGCAATCGACCGTTGCAGCAATATCTGGGCTTCACTGCCGGGCGCTGGTTATGGTCAGTTCGAGCATAAGGCTGACAGCCTGATTGCAAAATTCAAAGAAGCGGGCGGAACGGTCAGAGAGATTGATGTATGAGCAGAGTCACCGCGATTATCTCCGCTCTGGTTATCTGCATCATCGTCTGCCTGTCATGGGCTGTTAATCATTACCGTGATAACGCCATTACCTACAAAGCCCAGCGCGACAAAAATGCCAGAGAACTGAAGCTGGCGAACGCGGCAATTACTGACATGCAGATGCGTCAGCGTGATGTTGCTGCGCTCGATGCAAAATACACGAAGGAGTTAGCTGATGCTAAAGCTGAAAATGATGCTCTGCGTGATGATGTTGCCGCTGGTCGTCGTCGGTTGCACATCAAAGCAGTCTGTCAGTCAGTGCGTGAAGCCACCACCGCCTCCGGCGTGGATAATGCAGCCTCCCCCCGACTGGCAGACACCGCTGAACGGGATTATTTCACCCTCAGAGAGAGGCTGATCACTATGCAAAAACAACTGGAAGGAACCCAGAAGTATATTAATGAGCAGTGCAGATAGAGTTGCCCATATCGATGGGCAACTCATGCAATTATTGTGAGCAATACACACGCGCTTCCAGCGGAGTATAAATGCCTAAAGTAATAAAACCGAGCAATCCATTTACGAATGTTTGCTGGGTTTCTGTTTTAACAACATTTTCTGCGCCGCCACAAATTTTGGCTGCATCGACAGTTTTCTTCTGCCCAATTCCAGAAACGAAGAAATGATGGGTGATGGTTTCCTTTGGTGCTACTGCTGCCGGTTTGTTTTGAACAGTAAACGTCTGTTGAGCACATCCTGTAATAAGCAGGGCCAGCGCAGTAGCGAGTAGCATTTTTTTCATGGTGTTATTCCCGATGCTTTTTGAAGTTCGCAGAATCGTATGTGTAGAAAATTAAACAAACCCTAAACAATGAGTTGAAATTTCATATTGTTAATATTTATTAATGTATGTCAGGTGCGATGAATCGTCATTGTATTCCCGGATTAACTATGTCCACAGCCCTGACGGGGAACTTCTCTGCGGGAGTGTCCGGGAATAATTAAAACGATGCACACAGGGTTTAGCGCGTACACGTATTGCATTATGCCAACGCCCCGGTGCTGACACGGAAGAAACCGGACGTTATGATTTAGCGTGGAAAGATTTGTGTAGTGTTCTGAATGCTCTCAGTAAATAGTAATGAATTATCAAAGGTATAGTAATATCTTTTATGTTCATGGATATTTGTAACCCATCGGAAAACTCCTGCTTTAGCAAGATTTTCCCTGTATTGCTGAAATGTGATTTCTCTTGATTTCAACCTATCATAGGACGTTTCTATAAGATGCGTGTTTCTTGAGAATTTAACATTTACAACCTTTTTAAGTCCTTTTATTAACACGGTGTTATCGTTTTCTAACACGATGTGAATATTATCTGTGGCTAGATAGTAAATATAATGTGAGACGTTGTGACGTTTTAGTTCAGAATAAAACAATTCACAGTCTAAATCTTTTCGCACTTGATCGAATATTTCTTTAAAAATGGCAACCTGAGCCATTGGTAAAACCTTCCATGTGATACGAGGGCGCGTAGTTTGCATTATCGTTTTTATCGTTTCAATCTGGTCTGACCTCCTTGTGTTTTGTTGATGATTTATGTCAAATATTAGGAATGTTTTCACTTAATAGTATTGGTTGCGTAACAAAGTGCGGTCCTGCTGGCATTCTGGAGGGAAATACAACCGACAGATGTATGTAAGGCCAACGTGCTCAAATCTTCATACAGAAAGATTTGAAGTAATATTTTAACCGCTAGATGAAGAGCAAGCGCATGGAGCGACAAAATGAATAAAGAACAATCTGCTGATGATCCCTCCGTGGATCTGATTCGTGTAAAAAATATGCTTAATAGCACCATTTCTATGAGTTACCCTGATGTTGTAATTGCATGTATAGAACATAAGGTGTCTCTGGAAGCATTCAGAGCAATTGAGGCAGCGTTGGTGAAGCACGATAATAATATGAAGGATTATTCCCTGGTGGTTGACTGATCACCATAACTGCTAATCATTCAAACTATTTAGTCTGTGACAGAGCCAACACGCAGTCTGTCACTGTCAGGAAAGTGGTAAAACTGCAACTCAATTACTGCAATGCCCTCGTAATTAAGTGAATTTACAATATCGTCCTGTTCGGAGGGAAGAACGCGGGATGTTCATTCTTCATCACTTTTAATTGATGTATATGCTCTCTTTTCTGACGTTAGTCTCCGACGGCAGGCTTCAATGACCCAGGCTGAGAAATTCCCGGACCCTTTTTGCTCAAGAGCGATGTTAATTTGTTCAATCATTTGGTTAGGAAAGCGGATGTTGCGGGTTGTTGTTCTGCGGGTTCTGTTCTTCGTTGACATGAGGTTGCCCCGTATTCAGTGTCGCTGATTTGTATTGTCTGAAGTTGTTTTTACGTTAAGTTGATGCAGATCAATTAATACGATACCTGCGTCATAATTGATTATTTGACGTGGTTTGATGGCCTCCACGCACGTTGTGATATGTAGATGATAATCATTATCACTTTACGGGTCCTTTCCGGTGAAAAAAAAGGTACCAAAAAAAACATCGTCGTGAGTAGTGAACCGTAAGC
# <<< EMBEDDED_DCS_FASTA_END <<<
