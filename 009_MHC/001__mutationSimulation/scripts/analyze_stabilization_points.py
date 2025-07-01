#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from pathlib import Path

# === Parameters ===
INPUT_FOLDER = Path("results")
SUMMARY_SUFFIX = "_summary_all.tsv"
STABILIZATION_THRESHOLD = 0.05  # 5% change between bins
SCORE_COLUMNS = [
    "sift4g_score_mean",
    "polyphen2_hdiv_score_mean",
    "cadd_phred_mean",
    "frac_syn",
    "frac_missense",
    "frac_stop"
]

# === Output ===
stab_rows = []
summary_files = list(INPUT_FOLDER.glob("SBS*/SBS*_summary_all.tsv"))

for file in summary_files:
    signature = file.parent.name
    df = pd.read_csv(file, sep="\t")

    # Compute mean per mutation count
    df_grouped = df.groupby("n").agg({col: ['mean', 'std'] for col in SCORE_COLUMNS})
    df_grouped.columns = ['_'.join(col).strip() for col in df_grouped.columns.values]
    df_grouped = df_grouped.reset_index()

    for metric in SCORE_COLUMNS:
        means = df_grouped[["n", f"{metric}_mean"]].dropna().values
        for i in range(len(means) - 1):
            n1, v1 = means[i]
            n2, v2 = means[i + 1]
            if v1 == 0:
                continue
            rel_change = abs(v2 - v1) / abs(v1)
            if rel_change < STABILIZATION_THRESHOLD:
                stab_rows.append({
                    "signature": signature,
                    "metric": metric,
                    "stabilization_n": int(n2),
                    "final_mean": round(v2, 4)
                })
                break  # take the first time the metric stabilizes

# Save table
stab_df = pd.DataFrame(stab_rows)
stab_df.to_csv("results/stabilization_SBS_summary_all.tsv", sep="\t", index=False)

# === Plots ===
for metric in SCORE_COLUMNS:
    plt.figure(figsize=(6, 4))
    sns.histplot(stab_df[stab_df["metric"] == metric]["stabilization_n"], bins=10, kde=False)
    plt.title(f"Stabilization Point for {metric}")
    plt.xlabel("Mutation count where score stabilizes")
    plt.ylabel("Number of SBS signatures")
    plt.tight_layout()
    plt.savefig(f"results/stabilization_SBS_hist_{metric}.pdf")
    plt.close()





