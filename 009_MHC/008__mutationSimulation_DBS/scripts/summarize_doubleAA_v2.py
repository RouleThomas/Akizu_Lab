#!/usr/bin/env python3
import argparse, re
from pathlib import Path
import pandas as pd

def single_double_counts(parq: Path) -> tuple[int,int,int]:
    df = pd.read_parquet(parq, columns=["ref_codon", "codon_pos", "codon_index"])

    # coding rows only
    coding = df["ref_codon"].fillna("").astype(str).str.len() == 3
    df = df.loc[coding].copy()
    if df.empty:
        return 0, 0, 0

    if "codon_pos" in df.columns and df["codon_pos"].notna().any():
        cp = pd.to_numeric(df["codon_pos"], errors="coerce")
        # expected 0/1/2
        double = (cp == 2)
    else:
        # legacy fallback ONLY (not ideal, but keeps script usable)
        ci = pd.to_numeric(df["codon_index"], errors="coerce")
        max_ci = int(ci.max())
        double = (ci == 3) if max_ci >= 3 else (ci == 2)

    n_total = int(len(df))
    n_double = int(double.sum())
    n_single = n_total - n_double
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
    out_path = n_dir / "AA_change_summary.txt"
    out[["replicate","n_singleAA","n_doubleAA","n_total"]].to_csv(out_path, sep="\t", index=False)
    print(f"✅ wrote {out_path} ({len(out)} reps)")

def main(root: Path):
    for sig_dir in sorted(root.glob("*")):
        for n_dir in sorted(sig_dir.glob("n_*")):
            process_n_folder(n_dir)

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Summarize single vs double AA DBS events (uses codon_pos).")
    ap.add_argument("--root", required=True, help="Root like results_DBS")
    args = ap.parse_args()
    main(Path(args.root))
