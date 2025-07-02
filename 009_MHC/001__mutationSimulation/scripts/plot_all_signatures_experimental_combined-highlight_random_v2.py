#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem

# Parameters
INPUT_DIR = Path("results")
RANDOM_FILE = INPUT_DIR / "random_v2" / "random_v2_summary_all.tsv"
OUTPUT_PDF = "results/combined_signature_experimental_summary_errorbars-highlight_random_v2.pdf"

METRICS = [
    ("sift4g_score_mean", "SIFT4G Score", "Score (Lower = More Damaging)"),
    ("polyphen2_hdiv_score_mean", "PolyPhen-2 HDIV Score", "Score (Higher = More Damaging)"),
    ("cadd_phred_mean", "CADD PHRED Score", "Score (Higher = More Damaging)"),
    ("frac_syn", "Fraction Synonymous", "Fraction"),
    ("frac_missense", "Fraction Missense", "Fraction"),
    ("frac_stop", "Fraction Stop", "Fraction"),
]

# Load all *_summary_all.tsv files except those in SBS* folders
all_data = []
for file in sorted(INPUT_DIR.glob("*/*_summary_all.tsv")):
    folder_name = file.parent.name
    if folder_name.startswith("SBS"):
        continue  # Skip SBS folders
    df = pd.read_csv(file, sep="\t")
    df["signature"] = folder_name
    all_data.append(df)

df = pd.concat(all_data, ignore_index=True)

# Load random data
random_df = pd.read_csv(RANDOM_FILE, sep="\t")
random_df["signature"] = "random_flat96"

# Plotting
with PdfPages(OUTPUT_PDF) as pdf:
    for metric, title, ylabel in METRICS:
        plt.figure(figsize=(10, 10))

        # Plot all SBS signatures in black
        for sig in sorted(df["signature"].unique()):
            sub = df[df["signature"] == sig]
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
                     markersize=4, capsize=3, label="random_flat96")
        if x and y:
            plt.text(x[-1] + 1000, y[-1], "random_flat96", fontsize=8,
                     color="red", va='center', fontweight="bold")

        plt.title(title)
        plt.xlabel("Number of Mutations")
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        pdf.savefig()
        plt.close()

print(f"✅ Saved to {OUTPUT_PDF}")
