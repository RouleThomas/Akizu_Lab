#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem

# Parameters
INPUT_DIR = Path("results_contexts_v3")
OUTPUT_PDF = "results_contexts_v3/combined_signature_summary_errorbars_v4.pdf"
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
for file in sorted(INPUT_DIR.glob("*/*_summary_all_v4.tsv")):
    sig = file.parent.name
    df = pd.read_csv(file, sep="\t")
    df["signature"] = sig
    all_data.append(df)
df = pd.concat(all_data, ignore_index=True)

def add_threshold_bands(ax, metric):
    """
    Draw colored threshold bands and dashed reference lines
    for SIFT4G, PolyPhen2, and CADD.
    """
    cfg = {
        "sift4g_score_mean": {
            "lines": [0.05],
            "bands": [
                (0, 0.05, "Damaging", "#ff9999"),   # light red
                (0.05, None, "Benign", "#b3ffb3"),  # light green
            ],
        },
        "polyphen2_hdiv_score_mean": {
            "lines": [0.15, 0.85],
            "bands": [
                (0, 0.15, "Benign", "#b3ffb3"),                # green
                (0.15, 0.85, "Possibly damaging", "#b3d9ff"),  # blue
                (0.85, None, "Damaging", "#ff9999"),           # red
            ],
        },
        "cadd_phred_mean": {
            "lines": [2, 10, 20],
            "bands": [
                (0, 2, "Benign", "#b3ffb3"),             # green
                (2, 10, "Uncertain", "#b3d9ff"),         # blue
                (10, 20, "Likely damaging", "#ffcc99"),  # orange-red
                (20, None, "Damaging", "#ff9999"),       # red
            ],
        },
    }
    if metric not in cfg:
        return

    # Draw colored spans and dashed threshold lines
    ymin, ymax = ax.get_ylim()
    x_left, x_right = ax.get_xlim()
    for y in cfg[metric]["lines"]:
        ax.axhline(y, linestyle="--", color="k", linewidth=1, alpha=0.8)

    for low, high, label, color in cfg[metric]["bands"]:
        hi = ymax if high is None else high
        ax.axhspan(low, hi, facecolor=color, alpha=0.25, zorder=0)
        # label placement
        y_mid = (low + hi) / 2.0
        ax.text(x_right + 0.02*(x_right - x_left), y_mid, label,
                va="center", ha="left", fontsize=8, color="black")

    # expand right margin for labels
    ax.set_xlim(x_left, x_right + 0.10 * (x_right - x_left))


# Plotting
with PdfPages(OUTPUT_PDF) as pdf:
    for metric, title, ylabel in METRICS:
        fig, ax = plt.subplots(figsize=(10, 10))

        for sig in df["signature"].unique():
            sub = df[df["signature"] == sig]
            grouped = sub.groupby("n")[metric]
            x = sorted(grouped.groups.keys())
            y = [grouped.get_group(n).mean() for n in x]
            yerr = [sem(grouped.get_group(n)) for n in x]
            ax.errorbar(x, y, yerr=yerr, fmt='-o', label=sig,
                        linewidth=1, markersize=4, capsize=3)

            if x:
                ax.text(x[-1] + 0.03*(max(x)-min(x) if len(x)>1 else x[-1]),
                        y[-1], sig, fontsize=6, va='center')

        ax.set_title(title)
        ax.set_xlabel("Number of Mutations")
        ax.set_ylabel(ylabel)
        ax.grid(True)

        add_threshold_bands(ax, metric)

        fig.tight_layout(rect=[0, 0, 0.85, 1])
        ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5),
                  fontsize="x-small", ncol=2)

        pdf.savefig(fig)
        plt.close(fig)

print(f"✅ Saved to {OUTPUT_PDF}")
