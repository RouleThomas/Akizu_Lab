#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem

# -----------------------
# Params
# -----------------------
INPUT_DIR = Path("results_contexts")
OUTPUT_PDF = "results_contexts/combined_signature_summary_plots-highlight_ATCA.pdf"
HIGHLIGHT_SIG = "A[T>C]A"

METRICS = [
    ("sift4g_score_mean", "SIFT4G Score", "Score (Lower = More Damaging)"),
    ("polyphen2_hdiv_score_mean", "PolyPhen-2 HDIV Score", "Score (Higher = More Damaging)"),
    ("cadd_phred_mean", "CADD PHRED Score", "Score (Higher = More Damaging)"),
    ("frac_syn", "Fraction Synonymous", "Fraction"),
    ("frac_missense", "Fraction Missense", "Fraction"),
    ("frac_stop", "Fraction Stop", "Fraction"),
]

# -----------------------
# Load data
# -----------------------
all_data = []
for file in sorted(INPUT_DIR.glob("*/*_summary_all.tsv")):
    sig = file.parent.name
    df = pd.read_csv(file, sep="\t")
    df["signature"] = sig
    all_data.append(df)

if not all_data:
    raise SystemExit(f"No *_summary_all.tsv files found under {INPUT_DIR}")

sbs_df = pd.concat(all_data, ignore_index=True)

if HIGHLIGHT_SIG not in set(sbs_df["signature"]):
    raise SystemExit(f"Highlight signature '{HIGHLIGHT_SIG}' not found in data.")

highlight_df = sbs_df[sbs_df["signature"] == HIGHLIGHT_SIG].copy()

# -----------------------
# Plot
# -----------------------
with PdfPages(OUTPUT_PDF) as pdf:
    for metric, title, ylabel in METRICS:
        plt.figure(figsize=(10, 10))

        # Plot all *other* signatures in black
        for sig in sorted(sbs_df["signature"].unique()):
            if sig == HIGHLIGHT_SIG:
                continue
            sub = sbs_df[sbs_df["signature"] == sig]
            grouped = sub.groupby("n")[metric]
            if not grouped.groups:
                continue
            x = sorted(grouped.groups.keys())
            y = [grouped.get_group(n).mean() for n in x]
            yerr_vals = [sem(grouped.get_group(n)) for n in x]
            plt.errorbar(
                x, y, yerr=yerr_vals,
                fmt='-o', color='black', alpha=0.7,
                linewidth=1, markersize=3, capsize=3
            )
            if x and y:
                plt.text(x[-1] + 1000, y[-1], sig,
                         fontsize=6, color="black", va='center')

        # Plot highlighted signature in red
        grouped_h = highlight_df.groupby("n")[metric]
        xh = sorted(grouped_h.groups.keys())
        yh = [grouped_h.get_group(n).mean() for n in xh]
        yherr_vals = [sem(grouped_h.get_group(n)) for n in xh]
        plt.errorbar(
            xh, yh, yerr=yherr_vals,
            fmt='-o', color='red', linewidth=2,
            markersize=4, capsize=3, label=HIGHLIGHT_SIG
        )
        if xh and yh:
            plt.text(xh[-1] + 1000, yh[-1], HIGHLIGHT_SIG,
                     fontsize=8, color="red", va='center', fontweight="bold")

        plt.title(title)
        plt.xlabel("Number of Mutations")
        plt.ylabel(ylabel)
        plt.grid(True)
        plt.tight_layout(rect=[0, 0, 0.85, 1])
        pdf.savefig()
        plt.close()

print(f"✅ Saved to {OUTPUT_PDF}")
