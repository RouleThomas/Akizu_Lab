#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem

# Parameters
INPUT_DIR = Path("results_experimental")
OUTPUT_PDF = "results_experimental/combined_signature_summary_errorbars.pdf"
METRICS = [
    ("sift4g_score_mean", "SIFT4G Score", "Score (Lower = More Damaging)"),
    ("polyphen2_hdiv_score_mean", "PolyPhen-2 HDIV Score", "Score (Higher = More Damaging)"),
    ("cadd_phred_mean", "CADD PHRED Score", "Score (Higher = More Damaging)"),
    ("frac_syn", "Fraction Synonymous", "Fraction"),
    ("frac_missense", "Fraction Missense", "Fraction"),
    ("frac_stop", "Fraction Stop", "Fraction"),
]

# Load data
all_data = []
for file in sorted(INPUT_DIR.glob("*/*_summary_all.tsv")):
    sig = file.parent.name
    df = pd.read_csv(file, sep="\t")
    df["signature"] = sig
    all_data.append(df)
df = pd.concat(all_data, ignore_index=True)

# Plotting
with PdfPages(OUTPUT_PDF) as pdf:


    for metric, title, ylabel in METRICS:
        plt.figure(figsize=(10, 10))  # taller figure
        for sig in df["signature"].unique():
            sub = df[df["signature"] == sig]
            grouped = sub.groupby("n")[metric]
            x = sorted(grouped.groups.keys())
            y = [grouped.get_group(n).mean() for n in x]
            yerr = [sem(grouped.get_group(n)) for n in x]
            plt.errorbar(x, y, yerr=yerr, fmt='-o', label=sig, linewidth=1, markersize=4, capsize=3)

            # Add label at last x point
            try:
                plt.text(x[-1] + 1000, y[-1], sig, fontsize=6, va='center')  # adjust offset if needed
            except IndexError:
                continue

        plt.title(title)
        plt.xlabel("Number of Mutations")
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.tight_layout(rect=[0, 0, 0.85, 1])  # allow room for legend

        # Add 2-column legend outside
        plt.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize="x-small", ncol=2)
        pdf.savefig()
        plt.close()

    
print(f"✅ Saved to {OUTPUT_PDF}")
