#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem

# Parameters
INPUT_DIR = Path("results_experimental")
RANDOM_FILE = Path("results") / "Flat" / "Flat_summary_all.tsv"
OUTPUT_PDF = "results_experimental/combined_signature_summary_plots-highlight_random.pdf"

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
sbs_df = pd.concat(all_data, ignore_index=True)

# Load random data
random_df = pd.read_csv(RANDOM_FILE, sep="\t")
random_df["signature"] = "Flat"

# Plotting
with PdfPages(OUTPUT_PDF) as pdf:
    for metric, title, ylabel in METRICS:
        plt.figure(figsize=(10, 10))

        # Plot all SBS signatures in black
        for sig in sorted(sbs_df["signature"].unique()):
            sub = sbs_df[sbs_df["signature"] == sig]
            grouped = sub.groupby("n")[metric]
            x = sorted(grouped.groups.keys())
            y = [grouped.get_group(n).mean() for n in x]
            yerr = [sem(grouped.get_group(n)) for n in x]
            plt.errorbar(x, y, yerr=yerr, fmt='-o', color='black', alpha=0.7,
                         linewidth=1, markersize=3, capsize=3)
            if x and y:
                plt.text(x[-1] + 1000, y[-1], sig, fontsize=6,
                         color="black", va='center')

        # Plot random_v2 in red
        grouped = random_df.groupby("n")[metric]
        x = sorted(grouped.groups.keys())
        y = [grouped.get_group(n).mean() for n in x]
        yerr = [sem(grouped.get_group(n)) for n in x]
        plt.errorbar(x, y, yerr=yerr, fmt='-o', color='red', linewidth=2,
                     markersize=4, capsize=3, label="Flat")
        if x and y:
            plt.text(x[-1] + 1000, y[-1], "Flat", fontsize=8,
                     color="red", va='center', fontweight="bold")

        plt.title(title)
        plt.xlabel("Number of Mutations")
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        pdf.savefig()
        plt.close()

print(f"✅ Saved to {OUTPUT_PDF}")
