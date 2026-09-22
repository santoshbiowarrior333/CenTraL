#!/usr/bin/env python3
"""
analyze_ddct.py - ddCt analysis for the CenTraL Stage 1 qPCR.

Input: a CSV of Ct values with the columns
    sample,condition,target,strand,rt,replicate,ct
where
    sample     biological sample id (e.g. rep1, rep2, rep3)
    condition  treatment/state (one of them is the calibrator)
    target     e.g. chr2, chr17, GUSB
    strand     forward / reverse / na (use na for reference genes)
    rt         +RT or -RT
    replicate  technical replicate number
    ct         Ct value; leave empty or write NA for no amplification

What it does, matching the protocol:
    1. averages technical replicates per (sample, condition, target, strand, rt)
    2. checks the +RT to -RT gap (warns when below --min-rt-gap, default 5)
    3. dCt   = Ct(target, +RT) - Ct(reference, +RT) within the same sample+condition
    4. ddCt  = dCt - mean dCt of the calibrator condition (per target+strand)
    5. fold  = 2^-ddCt
Writes <out>/per_sample.tsv, <out>/summary.tsv and, if matplotlib is
available, <out>/fold_changes.png. Statistics (t test / ANOVA on dCt) are
left to your stats software; the per-sample dCt table is the input for that.

Usage:
    python3 analyze_ddct.py example_ct.csv --reference GUSB \
        --calibrator untreated --out results_qpcr
"""
import argparse
import csv
import math
import os
import sys
from collections import defaultdict


def mean(xs):
    return sum(xs) / len(xs)


def sd(xs):
    if len(xs) < 2:
        return float("nan")
    m = mean(xs)
    return math.sqrt(sum((x - m) ** 2 for x in xs) / (len(xs) - 1))


def norm_rt(v):
    v = v.strip().lower().replace(" ", "")
    if v in ("+rt", "plusrt", "rt+", "plus"):
        return "+RT"
    if v in ("-rt", "minusrt", "rt-", "minus", "nort", "no-rt"):
        return "-RT"
    sys.exit(f"ERROR: rt value not understood: '{v}' (use +RT or -RT)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("csv_file")
    ap.add_argument("--reference", default="GUSB",
                    help="reference target for dCt (default GUSB)")
    ap.add_argument("--calibrator", required=True,
                    help="condition used as the ddCt baseline")
    ap.add_argument("--min-rt-gap", type=float, default=5.0,
                    help="minimum acceptable +RT to -RT Ct gap (default 5)")
    ap.add_argument("--out", default="results_qpcr")
    args = ap.parse_args()

    rows = []
    with open(args.csv_file, newline="") as fh:
        reader = csv.DictReader(fh)
        need = {"sample", "condition", "target", "strand", "rt", "replicate", "ct"}
        missing = need - set(h.strip() for h in (reader.fieldnames or []))
        if missing:
            sys.exit(f"ERROR: missing columns: {sorted(missing)}")
        for r in reader:
            ct_raw = (r["ct"] or "").strip()
            ct = None if ct_raw.upper() in ("", "NA", "NAN", "UNDETERMINED") else float(ct_raw)
            rows.append({
                "sample": r["sample"].strip(),
                "condition": r["condition"].strip(),
                "target": r["target"].strip(),
                "strand": (r["strand"] or "na").strip() or "na",
                "rt": norm_rt(r["rt"]),
                "ct": ct,
            })
    if not rows:
        sys.exit("ERROR: no data rows")

    # 1. technical-replicate means
    tech = defaultdict(list)
    for r in rows:
        if r["ct"] is not None:
            tech[(r["sample"], r["condition"], r["target"], r["strand"], r["rt"])].append(r["ct"])
    ctmean = {k: mean(v) for k, v in tech.items()}

    # 2. RT gap check
    warnings = []
    for (smp, cond, tgt, strand, rt), m in sorted(ctmean.items()):
        if rt != "+RT":
            continue
        minus = ctmean.get((smp, cond, tgt, strand, "-RT"))
        if minus is None:
            continue
        gap = minus - m
        if gap < args.min_rt_gap:
            warnings.append(f"WARNING: +RT/-RT gap {gap:.1f} < {args.min_rt_gap} "
                            f"for {smp} {cond} {tgt} {strand} "
                            f"(possible gDNA or primer carryover)")

    # 3. per-sample dCt
    per_sample = []          # dicts: sample, condition, target, strand, ct, ref_ct, dct
    for (smp, cond, tgt, strand, rt), m in sorted(ctmean.items()):
        if rt != "+RT" or tgt == args.reference:
            continue
        refs = [v for (s2, c2, t2, _st, r2), v in ctmean.items()
                if s2 == smp and c2 == cond and t2 == args.reference and r2 == "+RT"]
        if not refs:
            warnings.append(f"WARNING: no {args.reference} +RT value for {smp} {cond}; "
                            f"skipping {tgt} {strand}")
            continue
        ref_ct = mean(refs)
        per_sample.append({"sample": smp, "condition": cond, "target": tgt,
                           "strand": strand, "ct": m, "ref_ct": ref_ct,
                           "dct": m - ref_ct})

    if not per_sample:
        sys.exit("ERROR: no dCt values could be computed")

    # 4-5. ddCt vs calibrator, per (target, strand)
    by_group = defaultdict(list)
    for p in per_sample:
        by_group[(p["target"], p["strand"], p["condition"])].append(p["dct"])
    calib_mean = {}
    for (tgt, strand, cond), dcts in by_group.items():
        if cond == args.calibrator:
            calib_mean[(tgt, strand)] = mean(dcts)
    summary = []
    for (tgt, strand, cond), dcts in sorted(by_group.items()):
        base = calib_mean.get((tgt, strand))
        if base is None:
            warnings.append(f"WARNING: no calibrator ({args.calibrator}) data for "
                            f"{tgt} {strand}; ddCt not computed for {cond}")
            continue
        ddct = mean(dcts) - base
        summary.append({"target": tgt, "strand": strand, "condition": cond,
                        "n": len(dcts), "mean_dct": mean(dcts), "sd_dct": sd(dcts),
                        "ddct": ddct, "fold_2^-ddct": 2 ** (-ddct)})
    for p in per_sample:                       # per-sample fold for plotting
        base = calib_mean.get((p["target"], p["strand"]))
        p["ddct"] = p["dct"] - base if base is not None else float("nan")
        p["fold"] = 2 ** (-p["ddct"]) if base is not None else float("nan")

    os.makedirs(args.out, exist_ok=True)
    with open(os.path.join(args.out, "per_sample.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["sample", "condition", "target", "strand", "mean_ct",
                    f"mean_ct_{args.reference}", "dct", "ddct", "fold"])
        for p in per_sample:
            w.writerow([p["sample"], p["condition"], p["target"], p["strand"],
                        f"{p['ct']:.3f}", f"{p['ref_ct']:.3f}", f"{p['dct']:.3f}",
                        f"{p['ddct']:.3f}", f"{p['fold']:.3f}"])
    with open(os.path.join(args.out, "summary.tsv"), "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["target", "strand", "condition", "n", "mean_dct", "sd_dct",
                    "ddct", "fold_2^-ddct"])
        for s in summary:
            w.writerow([s["target"], s["strand"], s["condition"], s["n"],
                        f"{s['mean_dct']:.3f}", f"{s['sd_dct']:.3f}",
                        f"{s['ddct']:.3f}", f"{s['fold_2^-ddct']:.3f}"])

    for wmsg in warnings:
        print(wmsg)
    print(f"wrote {args.out}/per_sample.tsv and {args.out}/summary.tsv")

    # optional plot
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("matplotlib not available; skipping plot")
        return
    groups = sorted({(s["target"], s["strand"]) for s in summary})
    conds = sorted({s["condition"] for s in summary},
                   key=lambda c: (c != args.calibrator, c))
    labels, heights, dots = [], [], []
    for tgt, strand in groups:
        for cond in conds:
            hit = [s for s in summary if s["target"] == tgt
                   and s["strand"] == strand and s["condition"] == cond]
            if not hit:
                continue
            labels.append(f"{tgt} {strand}\n{cond}")
            heights.append(hit[0]["fold_2^-ddct"])
            dots.append([p["fold"] for p in per_sample
                         if p["target"] == tgt and p["strand"] == strand
                         and p["condition"] == cond])
    fig, ax = plt.subplots(figsize=(max(4, 1.1 * len(labels)), 4))
    x = range(len(labels))
    ax.bar(x, heights, color="#c9d6e8", edgecolor="#3b6fc4", width=0.6)
    for xi, ys in zip(x, dots):
        ax.plot([xi] * len(ys), ys, "o", color="#333F4B", markersize=4, alpha=0.8)
    ax.axhline(1, linestyle="--", color="grey", linewidth=0.8)
    ax.set_xticks(list(x))
    ax.set_xticklabels(labels, fontsize=7)
    ax.set_ylabel("relative level (2^-ddCt)", fontsize=9)
    ax.set_title(f"CenTraL qPCR, calibrator: {args.calibrator}, "
                 f"reference: {args.reference}", fontsize=9)
    for spine in ("top", "right"):
        ax.spines[spine].set_visible(False)
    fig.tight_layout()
    fig.savefig(os.path.join(args.out, "fold_changes.png"), dpi=200)
    print(f"wrote {args.out}/fold_changes.png")


if __name__ == "__main__":
    main()
