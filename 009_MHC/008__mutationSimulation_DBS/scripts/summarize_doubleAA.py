#!/usr/bin/env python3

import argparse, re
from pathlib import Path
import pandas as pd

def single_double_counts(parq: Path) -> tuple[int,int,int]:
    """Return (n_total, n_singleAA, n_doubleAA) for one replicate parquet."""
    df = pd.read_parquet(parq, columns=["codon_index", "ref_codon"])
    # keep rows that look coding (have a 3-mer ref_codon)
    coding = df["ref_codon"].fillna("").astype(str).str.len() == 3
    df = df.loc[coding].copy()
    if df.empty:
        return 0, 0, 0

    # robust handling of 0-based vs 1-based indexing
    ci = pd.to_numeric(df["codon_index"], errors="coerce")
    max_ci = int(ci.max())
    if max_ci >= 3:
        # 1-based (1,2,3)
        double = (ci == 3)
    else:
        # 0-based (0,1,2)
        double = (ci == 2)

    n_total = int(len(ci))
    n_double = int(double.sum())
    n_single = int(n_total - n_double)
    return n_total, n_single, n_double

def process_n_folder(n_dir: Path):
    rows = []
    for p in sorted(n_dir.glob("rep_*.sim.parquet")):
        m = re.search(r"rep_(\d+)", p.stem)
        rep = f"rep{m.group(1).zfill(2)}" if m else p.stem
        n_total, n_single, n_double = single_double_counts(p)
        rows.append({"replicate": rep, "n_singleAA": n_single, "n_doubleAA": n_double, "n_total": n_total})

    if not rows:
        print(f"⚠️ No replicate Parquets in {n_dir}")
        return

    out = pd.DataFrame(rows).sort_values("replicate")
    # drop n_total if you only want the two columns; keeping it can help QC
    out_path = n_dir / "AA_change_summary.txt"
    out[["replicate","n_singleAA","n_doubleAA"]].to_csv(out_path, sep="\t", index=False)
    print(f"✅ wrote {out_path} ({len(out)} reps)")

def main(root: Path):
    # Walk results_*/*/n_*/
    for dbs_dir in sorted(root.glob("*")):
        for n_dir in sorted(dbs_dir.glob("n_*")):
            process_n_folder(n_dir)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Summarize single vs double AA events per replicate using codon_index.")
    ap.add_argument("--root", required=True, help="Root like results")
    args = ap.parse_args()
    main(Path(args.root))
