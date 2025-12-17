#!/usr/bin/env python3
from __future__ import annotations
import argparse, gzip
from pathlib import Path
import pandas as pd

def main():
    ap = argparse.ArgumentParser(description="Build a compact AA-level index from dbNSFP (for fast AA-change matching).")
    ap.add_argument("--dbnsfp", required=True, help="dbNSFP bgzip/gz file (e.g. ref/dbNSFP5.2a_grch38.gz)")
    ap.add_argument("--out", required=True, help="Output parquet (e.g. meta/dbnsfp_aa_index.parquet)")
    args = ap.parse_args()

    db = Path(args.dbnsfp)
    out = Path(args.out)
    out.parent.mkdir(parents=True, exist_ok=True)

    # Read header to get column indices by name (robust to version changes)
    with gzip.open(db, "rt") as fh:
        header = fh.readline().rstrip("\n").split("\t")

    name2idx = {c: i for i, c in enumerate(header)}

    # Required columns (names exactly as in dbNSFP header)
    required = [
        "#chr", "pos(1-based)", "ref", "alt",
        "aaref", "aaalt", "aapos",
        "Ensembl_geneid", "Ensembl_transcriptid",
        "SIFT4G_score", "Polyphen2_HDIV_score", "CADD_phred",
    ]
    missing = [c for c in required if c not in name2idx]
    if missing:
        raise SystemExit(f"❌ Missing columns in dbNSFP header: {missing}")

    use_idx = [name2idx[c] for c in required]

    rows = []
    with gzip.open(db, "rt") as fh:
        _ = fh.readline()  # skip header
        for line in fh:
            f = line.rstrip("\n").split("\t")
            rows.append([f[i] for i in use_idx])

    df = pd.DataFrame(rows, columns=required)

    # Normalize types
    df["pos(1-based)"] = pd.to_numeric(df["pos(1-based)"], errors="coerce")
    df["aapos"] = pd.to_numeric(df["aapos"], errors="coerce")

    for col in ["SIFT4G_score", "Polyphen2_HDIV_score", "CADD_phred"]:
        df[col] = pd.to_numeric(df[col], errors="coerce")

    df.to_parquet(out, index=False)
    print(f"✅ Wrote {out} with {len(df):,} rows")

if __name__ == "__main__":
    main()
