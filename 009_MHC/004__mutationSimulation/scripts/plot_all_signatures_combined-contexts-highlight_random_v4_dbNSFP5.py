#!/usr/bin/env python3

import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from pathlib import Path
from matplotlib.backends.backend_pdf import PdfPages
from scipy.stats import sem
import warnings

# --------------------- parameters ---------------------
INPUT_DIR = Path("results_contexts_v3")
# Try DBNSFP5 filename first, then plain v4 as fallback
FLAT_DBNSFP = Path("results_v3") / "Flat" / "Flat_summary_all_v4_dbNSFP5.tsv"
OUTPUT_PDF  = "results_contexts_v3/combined_signature_summary_plots-highlight_random_v4_dbNSFP5.pdf"

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

# --------------------- helpers ---------------------
def _safe_mean(s):
    return pd.to_numeric(s, errors="coerce").mean()

def _safe_sem(s):
    v = pd.to_numeric(s, errors="coerce").dropna()
    return sem(v) if len(v) > 1 else 0.0

def _load_many(dirpath: Path) -> pd.DataFrame:
    patterns = ["*/*_summary_all_v4_dbNSFP5.tsv", "*/*_summary_all_v4.tsv"]
    files = []
    for pat in patterns:
        files.extend(sorted(dirpath.glob(pat)))
    if not files:
        raise SystemExit(f"No summary TSVs found under {dirpath}")
    parts = []
    for f in files:
        sig = f.parent.name
        df = pd.read_csv(f, sep="\t")
        df["signature"] = sig
        parts.append(df)
    return pd.concat(parts, ignore_index=True)

def _load_flat() -> pd.DataFrame:
    if FLAT_DBNSFP.exists():
        f = FLAT_DBNSFP
    elif FLAT_V4.exists():
        f = FLAT_V4
    else:
        raise SystemExit("Flat summary TSV not found (tried DBNSFP5 and v4).")
    df = pd.read_csv(f, sep="\t")
    df["signature"] = "Flat"
    return df

# --------------------- load data ---------------------
sbs_df  = _load_many(INPUT_DIR)
flat_df = _load_flat()

# Warn about missing columns; skip those plots
missing = [c for c, _ in CATEGORIES if c not in sbs_df.columns or c not in flat_df.columns]
for m in missing:
    warnings.warn(f"Column missing; skipping plot: {m}")

categories = [(c, t) for c, t in CATEGORIES if c in sbs_df.columns and c in flat_df.columns]
if not categories:
    raise SystemExit("None of the expected percentage columns were found in both SBS and Flat data.")

# --------------------- plotting ---------------------
with PdfPages(OUTPUT_PDF) as pdf:
    for col, title in categories:
        fig, ax = plt.subplots(figsize=(10, 10))

        # All SBS signatures in black (mean ± SEM across reps at each n)
        for sig in sorted(sbs_df["signature"].unique()):
            sub = sbs_df[sbs_df["signature"] == sig]
            grp = sub.groupby("n")[col]
            if grp.ngroups == 0:
                continue
            ns = sorted(grp.groups.keys())
            y  = [_safe_mean(grp.get_group(n)) for n in ns]
            yerr = [_safe_sem(grp.get_group(n)) for n in ns]
            ax.errorbar(ns, y, yerr=yerr, fmt='-o', color='black', alpha=0.7,
                        linewidth=1, markersize=3, capsize=3)
            if ns and y:
                span = (max(ns) - min(ns)) if len(ns) > 1 else (ns[-1] or 1)
                ax.text(ns[-1] + 0.03 * span, y[-1], sig, fontsize=6,
                        color="black", va='center')

        # Flat highlighted in red
        grpF = flat_df.groupby("n")[col]
        if grpF.ngroups:
            nsF = sorted(grpF.groups.keys())
            yF  = [_safe_mean(grpF.get_group(n)) for n in nsF]
            yF_err = [_safe_sem(grpF.get_group(n)) for n in nsF]
            ax.errorbar(nsF, yF, yerr=yF_err, fmt='-o', color='red', linewidth=2,
                        markersize=4, capsize=3, label="Flat")
            if nsF and yF:
                spanF = (max(nsF) - min(nsF)) if len(nsF) > 1 else (nsF[-1] or 1)
                ax.text(nsF[-1] + 0.03 * spanF, yF[-1], "Flat", fontsize=8,
                        color="red", va='center', fontweight="bold")

        ax.set_title(title)
        ax.set_xlabel("Number of mutations (n)")
        ax.set_ylabel("Percent of variants")
        ax.grid(True, alpha=0.3)

        # expand right margin for end labels
        x_min, x_max = ax.get_xlim()
        ax.set_xlim(x_min, x_max + 0.10 * (x_max - x_min if x_max > x_min else 1))

        # Legend optional; end labels already present
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
