#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem
import warnings

# -------- Parameters --------
INPUT_DIR = Path("results_v3")
OUTPUT_PDF = "results_v3/combined_signature_summary_errorbars_v4_dbNSFP5.pdf"

# One plot per category (percentages 0–100)
CATEGORIES = [
    # SIFT4G (2)
    ("sift4g_damaging_pct", "SIFT4G — Damaging (≤ 0.05)"),
    ("sift4g_benign_pct",   "SIFT4G — Benign (> 0.05)"),
    # PolyPhen-2 HDIV (3)
    ("polyphen2_damaging_pct",          "PolyPhen-2 (HDIV) — Damaging (≥ 0.85)"),
    ("polyphen2_possibly_damaging_pct", "PolyPhen-2 (HDIV) — Possibly damaging (0.15–<0.85)"),
    ("polyphen2_benign_pct",            "PolyPhen-2 (HDIV) — Benign (< 0.15)"),
    # CADD PHRED (4)
    ("cadd_damaging_pct",        "CADD — Damaging (≥ 20)"),
    ("cadd_likely_damaging_pct", "CADD — Likely damaging (10–<20)"),
    ("cadd_uncertain_pct",       "CADD — Uncertain (2–<10)"),
    ("cadd_benign_pct",          "CADD — Benign (< 2)"),
]

# -------- Helpers --------
def _safe_mean(s):
    return pd.to_numeric(s, errors="coerce").mean()

def _safe_sem(s):
    v = pd.to_numeric(s, errors="coerce").dropna()
    return sem(v) if len(v) > 1 else 0.0

# -------- Load data --------
all_data = []

# Accept both file name variants
patterns = ["*/*_summary_all_v4_dbNSFP5.tsv"]
files = []
for pat in patterns:
    files.extend(sorted(INPUT_DIR.glob(pat)))

if not files:
    raise SystemExit(f"No summary TSVs found under {INPUT_DIR}")

for file in files:
    sig = file.parent.name
    df = pd.read_csv(file, sep="\t")
    df["signature"] = sig
    all_data.append(df)

df = pd.concat(all_data, ignore_index=True)

# Warn for missing columns (will skip those plots)
missing = [col for col, _ in CATEGORIES if col not in df.columns]
for m in missing:
    warnings.warn(f"Column missing; skipping plot: {m}")

categories = [(c, t) for c, t in CATEGORIES if c in df.columns]
if not categories:
    raise SystemExit("None of the expected percentage columns were found. Check your summarize script output.")

# -------- Plotting --------
with PdfPages(OUTPUT_PDF) as pdf:
    for col, title in categories:
        fig, ax = plt.subplots(figsize=(10, 10))

        for sig in sorted(df["signature"].unique()):
            sub = df[df["signature"] == sig]
            grouped = sub.groupby("n")[col]
            if grouped.ngroups == 0:
                continue
            ns = sorted(grouped.groups.keys())
            y_means = [_safe_mean(grouped.get_group(n)) for n in ns]
            y_errs  = [_safe_sem(grouped.get_group(n))  for n in ns]

            ax.errorbar(ns, y_means, yerr=y_errs, fmt='-o',
                        linewidth=1, markersize=4, capsize=3, label=sig)

            # Label last point for readability
            if ns:
                x_span = (max(ns) - min(ns)) if len(ns) > 1 else (ns[-1] or 1)
                ax.text(ns[-1] + 0.03 * x_span,
                        y_means[-1],
                        sig, fontsize=7, va='center')

        ax.set_title(title)
        ax.set_xlabel("Number of mutations (n)")
        ax.set_ylabel("Percent of variants")
        ax.grid(True, alpha=0.3)

        # Room on the right for end labels
        x_min, x_max = ax.get_xlim()
        ax.set_xlim(x_min, x_max + 0.10 * (x_max - x_min if x_max > x_min else 1))

        # Legend (optional; end labels already present)
        ax.legend(
            loc="center left",
            bbox_to_anchor=(1.02, 0.5),
            fontsize="x-small",
            ncol=2,           # ← two columns here
            columnspacing=0.8,
            handletextpad=0.4,
            borderaxespad=0.5
        )

        fig.tight_layout(rect=[0, 0, 0.88, 1])
        pdf.savefig(fig)
        plt.close(fig)

print(f"✅ Saved to {OUTPUT_PDF}")
