#!/usr/bin/env python3
"""Plot per-barcode read counts before and after DCS spike-in normalization.

Reads dcs_counts.tsv (output of count_dcs_spikein.sh) and writes a two-panel
PNG + PDF:
  - Top:    raw DCS spike-in counts per barcode (the basis of normalization)
  - Bottom: target-mapped reads, raw vs raw × scale_factor

Use this to sanity-check whether spike-in normalization actually equalizes
samples — bars in the bottom panel should be much more uniform after scaling
than before. If not, your input loading / sequencing was wildly uneven and
the spike-in is doing real work.

Usage:
    ./plot_normalization_qc.py -c /path/to/dcs_counts.tsv
    ./plot_normalization_qc.py -c dcs_counts.tsv -o my_qc_plot

Requires: matplotlib, numpy  (both in the python_collection env).
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # no display needed
import matplotlib.pyplot as plt
import numpy as np


def read_tsv(path):
    """Return list of dicts from a tab-delimited file with a header row."""
    rows = []
    with open(path) as f:
        header = f.readline().rstrip("\n").split("\t")
        for line_no, line in enumerate(f, start=2):
            parts = line.rstrip("\n").split("\t")
            if not parts or not parts[0]:
                continue
            rows.append(dict(zip(header, parts)))
    return rows


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    ap.add_argument("-c", "--counts", required=True,
                    help="Path to dcs_counts.tsv from count_dcs_spikein.sh")
    ap.add_argument("-o", "--out", default=None,
                    help="Output prefix (default: same dir as TSV, name 'dcs_normalization_qc')")
    args = ap.parse_args()

    tsv = Path(args.counts)
    if not tsv.exists():
        sys.exit(f"TSV not found: {tsv}")

    out_prefix = Path(args.out) if args.out else tsv.parent / "dcs_normalization_qc"

    try:
        rows = read_tsv(tsv)
    except Exception as e:
        sys.exit(f"Error reading TSV: {e}")

    if not rows:
        sys.exit("No data rows in TSV.")

    barcodes = [r["barcode"] for r in rows]
    
    # IMPROVED: Add error handling for numeric conversions
    try:
        target_mapped = np.array([int(r["target_mapped"]) for r in rows], dtype=float)
        dcs_mapped = np.array([int(r["dcs_mapped"]) for r in rows], dtype=float)
        scale = np.array(
            [float(r["scale_factor"]) if r["scale_factor"] != "NA" else np.nan
             for r in rows]
        )
    except (ValueError, KeyError) as e:
        sys.exit(f"Error parsing numeric columns: {e}")
    
    normalized = target_mapped * scale

    # Console summary so the user gets numbers as well as a picture
    print(f"Read {len(barcodes)} barcode(s) from {tsv}")
    try:
        print(f"  raw target_mapped — min: {int(target_mapped.min()):,}   "
              f"max: {int(target_mapped.max()):,}   "
              f"CV: {target_mapped.std()/target_mapped.mean():.2f}")
        finite_norm = normalized[np.isfinite(normalized)]
        if finite_norm.size:
            print(f"  normalized        — min: {int(finite_norm.min()):,}   "
                  f"max: {int(finite_norm.max()):,}   "
                  f"CV: {finite_norm.std()/finite_norm.mean():.2f}")
    except Exception as e:
        print(f"  Warning: Could not compute statistics: {e}", file=sys.stderr)

    # Plot — two panels, shared x axis
    width_per_bar = 0.6
    fig_w = max(8.0, len(barcodes) * width_per_bar)
    fig, (ax_top, ax_bot) = plt.subplots(2, 1, figsize=(fig_w, 8.5), sharex=True)

    # ---- top: DCS counts ----
    bar_color = "#6c757d"
    ax_top.bar(barcodes, dcs_mapped, color=bar_color, edgecolor="black", linewidth=0.4)
    ax_top.set_ylabel("DCS spike-in reads")
    ax_top.set_title("Per-barcode read counts — before vs after DCS spike-in normalization",
                     fontsize=12, pad=10)
    ax_top.grid(axis="y", alpha=0.25)
    for i, v in enumerate(dcs_mapped):
        ax_top.text(i, v, f"{int(v):,}", ha="center", va="bottom", fontsize=8)

    # ---- bottom: raw vs normalized ----
    x = np.arange(len(barcodes))
    bw = 0.4
    ax_bot.bar(x - bw/2, target_mapped, bw,
               label="raw target_mapped",
               color="#3a86ff", edgecolor="black", linewidth=0.4)
    ax_bot.bar(x + bw/2, normalized, bw,
               label="× scale_factor (DCS-normalized)",
               color="#06d6a0", edgecolor="black", linewidth=0.4)

    ax_bot.set_xticks(x)
    ax_bot.set_xticklabels(barcodes, rotation=45, ha="right")
    ax_bot.set_ylabel("read count")
    ax_bot.legend(loc="upper right", framealpha=0.9)
    ax_bot.grid(axis="y", alpha=0.25)

    for i, (raw, nrm) in enumerate(zip(target_mapped, normalized)):
        ax_bot.text(i - bw/2, raw, f"{int(raw/1000):,}k",
                    ha="center", va="bottom", fontsize=7)
        if np.isfinite(nrm):
            ax_bot.text(i + bw/2, nrm, f"{int(nrm/1000):,}k",
                        ha="center", va="bottom", fontsize=7)

    plt.tight_layout()
    
    try:
        png = out_prefix.with_suffix(".png")
        pdf = out_prefix.with_suffix(".pdf")
        plt.savefig(png, dpi=150)
        plt.savefig(pdf)
        print(f"\nWrote: {png}")
        print(f"Wrote: {pdf}")
    except Exception as e:
        sys.exit(f"Error writing output files: {e}")


if __name__ == "__main__":
    main()
