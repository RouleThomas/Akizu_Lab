#!/usr/bin/env python3

import argparse
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input", required=True, help="Path to signature summary TSV")
    parser.add_argument("--output", required=True, help="Output PDF path")
    args = parser.parse_args()

    df = pd.read_csv(args.input, sep="\t")
    grouped = df.groupby("n")

    # Metrics to plot: (column, label, ylabel, direction)
    metrics = [
        ("sift4g_score_mean", "SIFT4G Score", "Score (Lower = More Damaging)"),
        ("polyphen2_hdiv_score_mean", "PolyPhen-2 HDIV Score", "Score (Higher = More Damaging)"),
        ("cadd_phred_mean", "CADD PHRED Score", "Score (Higher = More Damaging)"),
        ("frac_syn", "Fraction Synonymous", "Fraction"),
        ("frac_missense", "Fraction Missense", "Fraction"),
        ("frac_stop", "Fraction Stop", "Fraction"),
    ]

    with PdfPages(args.output) as pdf:
        fig, axes = plt.subplots(3, 2, figsize=(12, 10))
        axes = axes.flatten()

        for i, (col, title, ylabel) in enumerate(metrics):
            ax = axes[i]
            x = sorted(grouped.groups.keys())
            y = [grouped.get_group(n)[col].mean() for n in x]
            yerr = [sem(grouped.get_group(n)[col]) for n in x]

            ax.errorbar(x, y, yerr=yerr, fmt='o-', capsize=4)
            ax.set_title(title)
            ax.set_xlabel("Number of Mutations")
            ax.set_ylabel(ylabel)
            ax.grid(True)

        plt.tight_layout()
        pdf.savefig(fig)
        plt.close()
        print(f"✅ Saved: {args.output}")

if __name__ == "__main__":
    main()
