import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path

# Define parameters
INPUT_DIR = Path("results")
METRICS = [
    ("sift4g_score_mean", "SIFT4G Score"),
    ("polyphen2_hdiv_score_mean", "PolyPhen2 HDIV Score"),
    ("cadd_phred_mean", "CADD PHRED Score"),
    ("frac_syn", "Fraction Synonymous"),
    ("frac_missense", "Fraction Missense"),
    ("frac_stop", "Fraction Stop"),
]

# Loop through all subdirectories in results/ that do NOT start with "SBS"
for subdir in sorted(INPUT_DIR.iterdir()):
    if not subdir.is_dir():
        continue
    if subdir.name.startswith("SBS"):
        continue

    summary_files = list(subdir.glob("*_summary_all.tsv"))
    if not summary_files:
        continue

    for file in summary_files:
        sig = subdir.name
        df = pd.read_csv(file, sep="\t")

        for metric, title in METRICS:
            if metric not in df.columns:
                print(f"⚠️ {metric} not found in {sig}, skipping.")
                continue

            grouped = df.groupby("n")[metric].agg(["mean", "std"]).reset_index()
            grouped["rel_change"] = grouped["mean"].pct_change().abs()

            # Find stabilization: 3 consecutive bins with <5% change
            stable = grouped["rel_change"] < 0.05
            rolling_stable = stable.rolling(window=3).apply(lambda x: x.all(), raw=True).fillna(0)
            stab_idx = rolling_stable.gt(0).idxmax()

            if stab_idx == 0 and not rolling_stable.any():
                print(f"❌ No stabilization found for {sig} - {metric}")
                stab_n = None
            else:
                stab_n = grouped.loc[stab_idx, "n"]
                print(f"📌 {sig}: {metric} stabilizes at ~{stab_n} mutations")

            # Plot
            plt.figure(figsize=(6, 4))
            plt.errorbar(grouped["n"], grouped["mean"], yerr=grouped["std"], fmt="-o", label=title)
            if stab_n:
                plt.axvline(stab_n, color="red", linestyle="--", label=f"Stabilized at {stab_n}")
            plt.title(f"{sig}: {title} Stabilization (3x <5%)")
            plt.xlabel("Number of Mutations")
            plt.ylabel(title)
            plt.legend()
            plt.tight_layout()

            # Save plot
            outfile = subdir / f"{sig}_{metric}_stabilization_check.pdf"
            plt.savefig(outfile)
            plt.close()


