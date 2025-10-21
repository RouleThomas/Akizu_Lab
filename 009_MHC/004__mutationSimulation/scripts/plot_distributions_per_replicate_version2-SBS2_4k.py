#!/usr/bin/env python3

import pandas as pd
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
from pathlib import Path
import numpy as np

# ------------------------- config -------------------------
BASE_DIR = Path("results_v3")

SCORE_COLUMNS = {
    "sift4g_score": "SIFT4G Score",
    "polyphen2_hdiv_score": "PolyPhen2 HDIV Score",
    "cadd_phred": "CADD PHRED Score",
    "frac_syn": "Fraction Synonymous",
    "frac_missense": "Fraction Missense",
    "frac_stop_gained": "Fraction Stop Gained",
}

# Category schemes (labels in display order + colors)
SCHEMES = {
    "sift4g_score": {
        "bins": [-np.inf, 0.05, np.inf],
        "labels": ["Damaging", "Benign"],
        "colors": ["#ff9999", "#b3ffb3"],  # red, green
    },
    "polyphen2_hdiv_score": {
        "bins": [-np.inf, 0.15, 0.85, np.inf],
        "labels": ["Benign", "Possibly damaging", "Damaging"],
        "colors": ["#b3ffb3", "#b3d9ff", "#ff9999"],  # green, blue, red
    },
    "cadd_phred": {
        "bins": [-np.inf, 2, 10, 20, np.inf],
        "labels": ["Benign", "Uncertain", "Likely damaging", "Damaging"],
        "colors": ["#b3ffb3", "#b3d9ff", "#ffcc99", "#ff9999"],  # green, blue, orange-red, red
    },
}

BINARY_COLORS = ["#d0d0d0", "#808080"]  # No, Yes
MISSING_COLOR = "#bfbfbf"               # gray for "Missing"

# ------------------------- helpers -------------------------
def consequence_to_binary(df: pd.DataFrame) -> pd.DataFrame:
    df["consequence"] = df["consequence"].astype(str)
    df["frac_syn"] = (df["consequence"] == "synonymous").astype("float")
    df["frac_missense"] = (df["consequence"] == "missense").astype("float")
    df["frac_stop_gained"] = (df["consequence"] == "stop").astype("float")
    # keep dtype float so any real NaN remains NaN if present
    return df

def categorized_counts(series: pd.Series, scheme_key: str):
    """
    Bucket numeric metric into qualitative categories and append a 'Missing' bar.
    Returns (labels, counts, colors).
    """
    scheme = SCHEMES[scheme_key]
    total = series.shape[0]
    observed = series.dropna()

    cats = pd.cut(observed, bins=scheme["bins"], labels=scheme["labels"], include_lowest=True)
    counts = cats.value_counts().reindex(scheme["labels"]).fillna(0).astype(int)

    missing = total - int(counts.sum())
    labels = list(scheme["labels"])
    counts = list(counts.values)
    colors = list(scheme["colors"])

    # Always add the Missing bar (0 if none) so totals are explicit
    labels.append("Missing")
    counts.append(missing)
    colors.append(MISSING_COLOR)

    return labels, counts, colors

def binary_counts(series: pd.Series):
    """
    Count Yes/No for a 0/1 series and append a 'Missing' bar.
    Returns (labels, counts, colors).
    """
    total = series.shape[0]
    observed = series.dropna().astype(int)

    # ensure both classes exist
    counts_series = observed.value_counts().reindex([0, 1]).fillna(0).astype(int)
    no_cnt = int(counts_series.loc[0]) if 0 in counts_series.index else 0
    yes_cnt = int(counts_series.loc[1]) if 1 in counts_series.index else 0
    missing = total - (no_cnt + yes_cnt)

    labels = ["No", "Yes", "Missing"]
    counts = [no_cnt, yes_cnt, missing]
    colors = BINARY_COLORS + [MISSING_COLOR]
    return labels, counts, colors

def add_bar_labels(ax, bars):
    for bar in bars:
        height = bar.get_height()
        ax.text(bar.get_x() + bar.get_width() / 2.0, height,
                f"{int(height)}", ha="center", va="bottom", fontsize=8)

# ------------------------- plotting -------------------------
# Iterate over the specific SBS folder (your example used SBS2/n_4000)
for sbs_folder in sorted(BASE_DIR.glob("SBS2/n_4000")):
    pdf_path = sbs_folder / f"{sbs_folder.parts[-2]}_replicate_score_distributions.pdf"
    with PdfPages(pdf_path) as pdf:
        # Find all replicates
        for parquet_file in sorted(sbs_folder.glob("rep_*.annot.parquet")):
            df = pd.read_parquet(parquet_file)
            df = consequence_to_binary(df)

            fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(8, 10))
            axs = axs.flatten()

            ax_idx = 0
            for col, title in SCORE_COLUMNS.items():
                ax = axs[ax_idx]
                ax_idx += 1

                if col not in df.columns:
                    ax.set_visible(False)
                    continue

                # Qualitative bar charts with "Missing"
                if col in SCHEMES:  # SIFT / PolyPhen / CADD
                    labels, counts, colors = categorized_counts(df[col], col)
                else:  # Binary consequence fractions -> Yes/No + Missing
                    labels, counts, colors = binary_counts(df[col])

                x = np.arange(len(labels))
                bars = ax.bar(x, counts, color=colors, edgecolor="black", linewidth=0.5)
                ax.set_xticks(x)
                ax.set_xticklabels(labels, rotation=20 if len(labels) > 3 else 0, ha="right")
                add_bar_labels(ax, bars)
                ax.set_ylabel("Count")
                ax.set_title(title)
                ax.spines["top"].set_visible(False)
                ax.spines["right"].set_visible(False)
                ax.grid(axis="y", linestyle=":", alpha=0.5)

            # Hide any unused subplots (if any)
            for j in range(ax_idx, len(axs)):
                axs[j].set_visible(False)

            replicate_name = parquet_file.stem
            total_mut = len(df)
            fig.suptitle(
                f"{sbs_folder.parts[-2]} - {replicate_name}: Mutation Score Categories (N={total_mut})",
                fontsize=12,
            )
            plt.tight_layout(rect=[0, 0.03, 1, 0.95])
            pdf.savefig(fig)
            plt.close(fig)

    print(f"✅ Finished plotting {sbs_folder.parts[-2]}")
