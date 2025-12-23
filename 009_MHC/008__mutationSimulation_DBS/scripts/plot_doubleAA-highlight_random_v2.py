#!/usr/bin/env python3
import pandas as pd
import matplotlib.pyplot as plt
from pathlib import Path
import re
import argparse

def _collect_from_sig_dir(sig_dir: Path, sig_name: str) -> list[pd.DataFrame]:
    recs = []
    for n_dir in sorted(sig_dir.glob("n_*")):
        m = re.search(r"n_(\d+)", n_dir.name)
        n_mut = int(m.group(1)) if m else None
        f = n_dir / "AA_change_summary.txt"
        if not f.exists():
            continue
        df = pd.read_csv(f, sep="\t")
        df["signature"] = sig_name
        df["n"] = n_mut
        df["prop_doubleAA"] = df["n_doubleAA"] / (df["n_singleAA"] + df["n_doubleAA"]).clip(lower=1)
        recs.append(df)
    return recs

def load_all_summaries(root: Path, flat_dir: Path | None) -> pd.DataFrame:
    records: list[pd.DataFrame] = []

    # * come from root only
    for sig_dir in sorted(root.glob("*")):
        if sig_dir.is_dir():
            records += _collect_from_sig_dir(sig_dir, sig_dir.name)

    # Flat comes from explicit path if provided
    if flat_dir is not None and flat_dir.is_dir():
        records += _collect_from_sig_dir(flat_dir, "Flat")

    return pd.concat(records, ignore_index=True) if records else pd.DataFrame()

def plot_prop(df: pd.DataFrame, outpath: Path, title: str | None = None):
    plt.figure(figsize=(9, 9))

    agg = (
        df.groupby(["signature", "n"])["prop_doubleAA"]
          .agg(["mean", "std"])
          .reset_index()
    )
    max_n = agg["n"].max()

    # Non-Flat (thin gray)
    others = agg[agg["signature"] != "Flat"]
    for sig, sub in others.groupby("signature"):
        plt.errorbar(sub["n"], sub["mean"]*100, yerr=sub["std"]*100,
                     color="0.35", alpha=0.8, linewidth=1.0, marker='o',
                     markersize=3, capsize=2, zorder=1)
        last_row = sub.loc[sub["n"].idxmax()]
        plt.text(max_n * 1.01, last_row["mean"]*100, sig, color="black",
                 fontsize=7, va='center', ha='left')

    # Flat highlighted in red (only if present)
    flat = agg[agg["signature"] == "Flat"]
    if not flat.empty:
        plt.errorbar(flat["n"], flat["mean"]*100, yerr=flat["std"]*100,
                     color="red", linewidth=2.5, marker='s', markersize=6,
                     capsize=3, label="Flat", zorder=3)
        last_flat = flat.loc[flat["n"].idxmax()]
        plt.text(max_n * 1.01, last_flat["mean"]*100, "Flat", color="red",
                 fontsize=8, fontweight="bold", va='center', ha='left')

    plt.xlabel("Number of simulated mutations (n)")
    plt.ylabel("Percent of double-AA mutations")
    plt.title(title or "Proportion of double-AA DBS events per signature")
    plt.xlim(0, max_n * 1.1)
    plt.tight_layout()
    plt.savefig(outpath, dpi=300, bbox_inches="tight")
    print(f"✅ Saved plot with labels to {outpath}")

def main():
    ap = argparse.ArgumentParser(
        description="Plot % double-AA vs n. DBS* loaded from --root; Flat loaded from --flat (independent)."
    )
    ap.add_argument("--root", required=True,
                    help="Root with DBS*/n_*/AA_change_summary.txt")
    ap.add_argument("--flat", default=None,
                    help="Path to Flat directory (…/Flat with n_*/AA_change_summary.txt). Omit to skip Flat.")
    ap.add_argument("--out", default="DBS_doubleAA_proportion_labeled.pdf", help="Output PDF")
    ap.add_argument("--title", default=None, help="Optional plot title")
    args = ap.parse_args()

    root = Path(args.root)
    flat_dir = Path(args.flat) if args.flat else None

    df = load_all_summaries(root, flat_dir)
    if df.empty:
        raise SystemExit("❌ No AA_change_summary.txt files found for DBS* and/or Flat (if provided).")

    plot_prop(df, Path(args.out), args.title)

if __name__ == "__main__":
    main()
