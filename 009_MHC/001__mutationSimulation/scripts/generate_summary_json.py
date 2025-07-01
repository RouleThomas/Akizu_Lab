import pandas as pd
import numpy as np
import argparse
import json
from pathlib import Path

def summarize_output(df):
    consequence_counts = df["consequence"].value_counts().to_dict()
    scores = {
        "sift4g_score": df["sift4g_score"].dropna().tolist(),
        "polyphen2_hdiv_score": df["polyphen2_hdiv_score"].dropna().tolist(),
        "cadd_phred": df["cadd_phred"].dropna().tolist(),
    }
    return {
        "n_mutations": len(df),
        "consequences": consequence_counts,
        "score_counts": {k: len(v) for k, v in scores.items()},
        "score_means": {k: round(np.mean(v), 3) if v else None for k, v in scores.items()},
    }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-dir", required=True, help="Path to the base folder (e.g., results/random)")
    args = parser.parse_args()

    base_dir = Path(args.input_dir)
    parquet_files = base_dir.glob("n_*/rep_*.annot.parquet")

    for pfile in parquet_files:
        json_path = pfile.with_suffix(".summary.json")
        df = pd.read_parquet(pfile)
        summary = summarize_output(df)
        with open(json_path, "w") as f:
            json.dump(summary, f, indent=2)
        print(f"✅ Created: {json_path}")

if __name__ == "__main__":
    main()
