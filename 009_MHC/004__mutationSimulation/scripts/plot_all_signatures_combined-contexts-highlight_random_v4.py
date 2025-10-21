#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem

# --------------------- thresholds helper ---------------------
def add_threshold_bands(ax, metric):
    """
    Draw colored bands + dashed reference lines for SIFT4G, PolyPhen2 (HDIV), and CADD.
    Colors: green=benign, blue=uncertain/possibly, red=damaging.
    """
    cfg = {
        "sift4g_score_mean": {
            "lines": [0.05],
            "bands": [
                (0, 0.05, "Damaging", "#ff9999"),    # red
                (0.05, None, "Benign", "#b3ffb3"),   # green
            ],
        },
        "polyphen2_hdiv_score_mean": {
            "lines": [0.15, 0.85],
            "bands": [
                (0, 0.15, "Benign", "#b3ffb3"),                 # green
                (0.15, 0.85, "Possibly damaging", "#b3d9ff"),   # blue
                (0.85, None, "Damaging", "#ff9999"),            # red
            ],
        },
        "cadd_phred_mean": {
            "lines": [2, 10, 20],
            "bands": [
                (0, 2, "Benign", "#b3ffb3"),              # green
                (2, 10, "Uncertain", "#b3d9ff"),          # blue
                (10, 20, "Likely damaging", "#ffcc99"),   # orange-red
                (20, None, "Damaging", "#ff9999"),        # red
            ],
        },
    }
    if metric not in cfg:
        return

    # dashed lines
    for y in cfg[metric]["lines"]:
        ax.axhline(y, linestyle="--", color="k", linewidth=1, alpha=0.85, zorder=1)

    # spans + right-hand labels
    ymin, ymax = ax.get_ylim()
    x_left, x_right = ax.get_xlim()
    for low, high, label, color in cfg[metric]["bands"]:
        hi = ymax if high is None else high
        ax.axhspan(low, hi, facecolor=color, alpha=0.25, zorder=0)
        y_mid = (low + hi) / 2.0
        ax.text(x_right + 0.02 * (x_right - x_left), y_mid, label,
                va="center", ha="left", fontsize=8, color="black")

    # add a little right margin for the labels
    ax.set_xlim(x_left, x_right + 0.10 * (x_right - x_left))

# --------------------- parameters ---------------------
INPUT_DIR = Path("results_contexts_v3")
RANDOM_FILE = Path("results_v3") / "Flat" / "Flat_summary_all_v4.tsv"
OUTPUT_PDF = "results_contexts_v3/combined_signature_summary_plots-highlight_random_v4.pdf"

METRICS = [
    ("sift4g_score_mean", "SIFT4G Score", "Score (Lower = More Damaging)"),
    ("polyphen2_hdiv_score_mean", "PolyPhen-2 HDIV Score", "Score (Higher = More Damaging)"),
    ("cadd_phred_mean", "CADD PHRED Score", "Score (Higher = More Damaging)"),
    ("frac_syn", "Fraction Synonymous", "Fraction"),
    ("frac_missense", "Fraction Missense", "Fraction"),
    ("frac_stop", "Fraction Stop", "Fraction"),
]

# --------------------- load data ---------------------
all_data = []
for file in sorted(INPUT_DIR.glob("*/*_summary_all_v4.tsv")):
    sig = file.parent.name
    df = pd.read_csv(file, sep="\t")
    df["signature"] = sig
    all_data.append(df)
sbs_df = pd.concat(all_data, ignore_index=True)

random_df = pd.read_csv(RANDOM_FILE, sep="\t")
random_df["signature"] = "Flat"

# --------------------- plotting ---------------------
with PdfPages(OUTPUT_PDF) as pdf:
    for metric, title, ylabel in METRICS:
        fig, ax = plt.subplots(figsize=(10, 10))

        # Plot all SBS signatures in black
        for sig in sorted(sbs_df["signature"].unique()):
            sub = sbs_df[sbs_df["signature"] == sig]
            grouped = sub.groupby("n")[metric]
            if not grouped.groups:
                continue
            x = sorted(grouped.groups.keys())
            y = [grouped.get_group(n).mean() for n in x]
            yerr = [sem(grouped.get_group(n)) for n in x]
            ax.errorbar(x, y, yerr=yerr, fmt='-o', color='black', alpha=0.7,
                        linewidth=1, markersize=3, capsize=3)
            if x and y:
                ax.text(x[-1] + 1000, y[-1], sig, fontsize=6,
                        color="black", va='center')

        # Plot Flat in red (highlight)
        grouped = random_df.groupby("n")[metric]
        if grouped.groups:
            x = sorted(grouped.groups.keys())
            y = [grouped.get_group(n).mean() for n in x]
            yerr = [sem(grouped.get_group(n)) for n in x]
            ax.errorbar(x, y, yerr=yerr, fmt='-o', color='red', linewidth=2,
                        markersize=4, capsize=3, label="Flat")
            if x and y:
                ax.text(x[-1] + 1000, y[-1], "Flat", fontsize=8,
                        color="red", va='center', fontweight="bold")

        ax.set_title(title)
        ax.set_xlabel("Number of Mutations")
        ax.set_ylabel(ylabel)
        ax.grid(True)

        # Add benign/possibly/damaging bands when applicable
        add_threshold_bands(ax, metric)

        fig.tight_layout(rect=[0, 0, 0.85, 1])  # keep room for right-side labels
        pdf.savefig(fig)
        plt.close(fig)

print(f"✅ Saved to {OUTPUT_PDF}")
