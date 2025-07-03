import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Input directory and metrics
INPUT_DIR = Path("results")
METRICS = {
    "sift4g_score_mean": "SIFT4G Score",
    "polyphen2_hdiv_score_mean": "PolyPhen2 HDIV Score",
    "cadd_phred_mean": "CADD PHRED Score",
    "frac_syn": "Fraction Synonymous",
    "frac_missense": "Fraction Missense",
    "frac_stop": "Fraction Stop",
}

# Store stabilization points
stabilization_points = {metric: {} for metric in METRICS}

# Extract stabilization points per SBS signature and metric
for sig_folder in sorted(INPUT_DIR.glob("SBS*")):
    sig_name = sig_folder.name
    summary_file = sig_folder / f"{sig_name}_summary_all.tsv"
    if not summary_file.exists():
        continue

    df = pd.read_csv(summary_file, sep="\t")
    for metric in METRICS:
        if metric not in df.columns:
            continue

        grouped = df.groupby("n")[metric].agg(["mean", "std"]).reset_index()
        grouped["rel_change"] = grouped["mean"].pct_change().abs()
        stable = grouped["rel_change"] < 0.05
        rolling_stable = stable.rolling(window=3).apply(lambda x: x.all(), raw=True).fillna(0)

        if rolling_stable.gt(0).any():
            stab_idx = rolling_stable.idxmax()
            stab_n = grouped.loc[stab_idx, "n"]
            stabilization_points[metric][sig_name] = stab_n

# Gather all possible stabilization n values across all metrics
all_stab_ns = sorted(set(
    stab_n
    for metric_dict in stabilization_points.values()
    for stab_n in metric_dict.values()
))

# 📊 Plot per metric with unified x-axis
for metric, stab_dict in stabilization_points.items():
    if not stab_dict:
        continue

    # Count and reindex with all possible values (fill missing with 0)
    stab_series = pd.Series(stab_dict).map(lambda x: x if x in all_stab_ns else None)
    value_counts = stab_series.value_counts().reindex(all_stab_ns, fill_value=0)

    # Plot
    plt.figure(figsize=(8, 5))
    bars = plt.bar(value_counts.index.astype(str), value_counts.values, color="navy", edgecolor="black")

    # Add count on top of each bar
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            plt.text(bar.get_x() + bar.get_width() / 2, height + 0.5, f"{int(height)}",
                     ha='center', va='bottom', fontsize=9)

    plt.title(f"Stabilization Point Distribution: {METRICS[metric]}")
    plt.xlabel("Number of Mutations at Stabilization")
    plt.ylabel("Number of SBS Signatures")
    plt.xticks(rotation=45)
    plt.tight_layout()

    plt.savefig(INPUT_DIR / f"stabilization_summary_SBS_{metric}.pdf")
    plt.close()

