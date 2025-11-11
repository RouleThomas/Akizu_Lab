#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem
import numpy as np
import warnings

def _safe_mean(series):
    return pd.to_numeric(series, errors="coerce").mean()

def _safe_sem(series):
    vals = pd.to_numeric(series, errors="coerce").dropna()
    return sem(vals) if len(vals) > 1 else 0.0

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to signature summary TSV")
    parser.add_argument("--output", required=True, help="Output PDF path")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")

    # One plot per category (percentages 0–100)
    metrics = [
        # SIFT4G
        ("sift4g_damaging_pct", "SIFT4G — Damaging (≤ 0.05)", "Percent of variants"),
        ("sift4g_benign_pct",   "SIFT4G — Benign (> 0.05)",   "Percent of variants"),
        # PolyPhen-2 (HDIV)
        ("polyphen2_damaging_pct",          "PolyPhen-2 (HDIV) — Damaging (≥ 0.85)",             "Percent of variants"),
        ("polyphen2_possibly_damaging_pct", "PolyPhen-2 (HDIV) — Possibly damaging (0.15–<0.85)", "Percent of variants"),
        ("polyphen2_benign_pct",            "PolyPhen-2 (HDIV) — Benign (< 0.15)",               "Percent of variants"),
        # CADD PHRED
        ("cadd_damaging_pct",        "CADD — Damaging (≥ 20)",          "Percent of variants"),
        ("cadd_likely_damaging_pct", "CADD — Likely damaging (10–<20)", "Percent of variants"),
        ("cadd_uncertain_pct",       "CADD — Uncertain (2–<10)",        "Percent of variants"),
        ("cadd_benign_pct",          "CADD — Benign (< 2)",             "Percent of variants"),
    ]

    # warn if columns are missing; they will be skipped
    missing = [col for col, *_ in metrics if col not in df.columns]
    for col in missing:
        warnings.warn(f"Column missing in input and will be skipped: {col}")

    metrics = [m for m in metrics if m[0] in df.columns]
    if not metrics:
        raise SystemExit("No expected metric columns found. Check your input TSV.")

    grouped = df.groupby("n")
    ns = sorted(grouped.groups.keys())

    # layout: make a square-ish grid
    nplots = len(metrics)
    ncols = 3
    nrows = int(np.ceil(nplots / ncols))

    with PdfPages(args.output) as pdf:
        fig, axes = plt.subplots(nrows, ncols, figsize=(4 * ncols + 2, 3.6 * nrows + 1.5))
        if nrows == 1 and ncols == 1:
            axes = np.array([[axes]])
        elif nrows == 1:
            axes = np.array([axes])
        axes = axes.flatten()

        for i, (col, title, ylabel) in enumerate(metrics):
            ax = axes[i]
            y_means = [_safe_mean(grouped.get_group(n)[col]) for n in ns]
            y_errs  = [_safe_sem(grouped.get_group(n)[col])  for n in ns]

            ax.errorbar(ns, y_means, yerr=y_errs, fmt='o-', capsize=4)
            ax.set_title(title, fontsize=11)
            ax.set_xlabel("Number of mutations (n)")
            ax.set_ylabel(ylabel)
            ax.set_ylim(0, 100)  # percentages
            ax.grid(True, alpha=0.3)

        # hide any unused axes
        for j in range(i + 1, len(axes)):
            axes[j].axis('off')

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close(fig)

    print(f"✅ Saved: {args.output}")

if __name__ == "__main__":
    main()
