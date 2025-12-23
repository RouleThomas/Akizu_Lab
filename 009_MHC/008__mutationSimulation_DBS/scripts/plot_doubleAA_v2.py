#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re
import argparse

def load_all_summaries(root: Path):
    records = []
    for dbs_dir in sorted(root.glob("*")):
        if not dbs_dir.is_dir():
            continue
        dbs = dbs_dir.name
        for n_dir in sorted(dbs_dir.glob("n_*")):
            m = re.search(r"n_(\d+)", n_dir.name)
            n_mut = int(m.group(1)) if m else None
            summary = n_dir / "AA_change_summary.txt"
            if summary.exists():
                df = pd.read_csv(summary, sep="\t")
                df["DBS"] = dbs
                df["n"] = n_mut
                df["prop_doubleAA"] = df["n_doubleAA"] / (df["n_singleAA"] + df["n_doubleAA"]).clip(lower=1)
                records.append(df)
    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()

def plot_prop(df: pd.DataFrame, outpath: Path):
    plt.figure(figsize=(9, 7))

    agg = (
        df.groupby(["DBS", "n"])["prop_doubleAA"]
          .agg(["mean", "std"])
          .reset_index()
    )

    max_n = agg["n"].max()

    for dbs, sub in agg.groupby("DBS"):
        plt.errorbar(sub["n"], sub["mean"]*100, yerr=sub["std"]*100,
                     marker='o', capsize=3, label=dbs, alpha=0.9)
        # Add label on right at max_n
        last = sub.loc[sub["n"].idxmax()]
        plt.text(max_n * 1.01, last["mean"]*100, dbs,
                 fontsize=7, color="black", va='center', ha='left')

    plt.xlabel("Number of simulated mutations (n)")
    plt.ylabel("Percent of double-AA mutations")
    plt.title("Proportion of double-AA DBS events per signature")
    plt.xlim(0, max_n * 1.1)  # space for labels on right
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"✅ Saved plot with DBS labels to {outpath}")

def main():
    parser = argparse.ArgumentParser(description="Plot % of double-AA DBS events per signature (with right-side labels).")
    parser.add_argument("--root", required=True, help="Root directory containing DBS*/n_*/AA_change_summary.txt")
    parser.add_argument("--out", default="DBS_doubleAA_proportion_labeled.pdf", help="Output plot path (PDF)")
    args = parser.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"❌ Root path not found: {root}")

    df = load_all_summaries(root)
    if df.empty:
        print(f"⚠️ No data found in {root}")
        return

    plot_prop(df, Path(args.out))

if __name__ == "__main__":
    main()
