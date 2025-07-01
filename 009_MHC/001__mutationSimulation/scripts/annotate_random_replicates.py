#!/usr/bin/env python3
"""
Annotate random flat96 simulations with damage scores and create JSON summaries.
"""

import argparse, json, sys
from pathlib import Path
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
import numpy as np
from annotate_damage2 import DBNSFP

def summarize_output(df):
    consequence_counts = df["consequence"].value_counts().to_dict()
    scores = {
        "sift4g_score": df.get("sift4g_score", pd.Series(dtype=float)).dropna().tolist(),
        "polyphen2_hdiv_score": df.get("polyphen2_hdiv_score", pd.Series(dtype=float)).dropna().tolist(),
        "cadd_phred": df.get("cadd_phred", pd.Series(dtype=float)).dropna().tolist(),
    }
    return {
        "n_mutations": len(df),
        "consequences": consequence_counts,
        "score_counts": {k: len(v) for k, v in scores.items()},
        "score_means": {k: round(np.mean(v), 3) if v else None for k, v in scores.items()},
    }

def annotate_file(parquet_path: Path, dbnsfp_path: Path):
    df = pq.read_table(parquet_path).to_pandas()

    db = DBNSFP(dbnsfp_path)
    annotations = {k: [] for k in [
        "polyphen2_hdiv_score", "polyphen2_hdiv_pred",
        "sift4g_score", "sift4g_pred", "cadd_phred"]}

    for chrom, pos, r_idx, a in zip(df.chr, df.pos, df.ref_base, df.alt_base):
        ref = "ACGT"[int(r_idx)]
        result = db.query(str(chrom), int(pos), ref, str(a))
        for key in annotations:
            annotations[key].append(result[key])

    for k, v in annotations.items():
        df[k] = v

    out_annot = parquet_path.with_name(parquet_path.stem + ".annot.parquet")
    out_json = parquet_path.with_name(parquet_path.stem + ".summary.json")

    pq.write_table(pa.Table.from_pandas(df), out_annot, compression="zstd")
    with open(out_json, "w") as f:
        json.dump(summarize_output(df), f, indent=2)

    print(f"✅ Annotated: {parquet_path.name}")

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, help="Folder with rep_*.parquet files")
    parser.add_argument("--dbnsfp", required=True, help="Path to dbNSFP file")
    args = parser.parse_args()

    input_dir = Path(args.input_dir)
    if not input_dir.exists():
        sys.exit(f"❌ ERROR: Input directory does not exist: {input_dir}")
    
    for file in sorted(input_dir.glob("rep_*.parquet")):
        if ".annot" in file.stem:
            continue  # Skip already-annotated files
        annotate_file(file, Path(args.dbnsfp))

if __name__ == "__main__":
    main()
