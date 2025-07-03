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

# Get all x values to ensure consistent x-axis
all_ns = sorted(set(
    val for stab_dict in stabilization_points.values() for val in stab_dict.values()
))

# Create 2-column, 3-row subplot grid
fig, axs = plt.subplots(nrows=3, ncols=2, figsize=(5, 8))
axs = axs.flatten()

for ax, (metric, title) in zip(axs, METRICS.items()):
    counts = pd.Series(stabilization_points[metric]).value_counts()
    counts = counts.reindex(all_ns, fill_value=0).sort_index()

    bars = ax.bar(counts.index.astype(str), counts.values, color="navy", edgecolor="black")

    # Add count on top of each bar
    for bar in bars:
        height = bar.get_height()
        if height > 0:
            ax.text(bar.get_x() + bar.get_width() / 2, height + 0.5, f"{int(height)}",
                    ha='center', va='bottom', fontsize=8)

    ax.set_title(f"Stabilization Point: {title}", fontsize=10)
    ax.set_ylabel("SBS Count", fontsize=9)
    ax.set_xlabel("Mutations at Stabilization", fontsize=9)
    ax.tick_params(axis='x', labelrotation=45)

# Remove unused subplots if fewer than total
for ax in axs[len(METRICS):]:
    ax.axis("off")

plt.tight_layout()
plt.savefig(INPUT_DIR / "stabilization_summary_SBS_all.pdf")
plt.close()


