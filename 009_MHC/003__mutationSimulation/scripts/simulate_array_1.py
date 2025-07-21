#!/usr/bin/env python3
"""
Wrapper script for a single SBS signature simulation + annotation run.
"""

import argparse, json, subprocess
from pathlib import Path
import pyarrow.parquet as pq


def summarize_output(df):
    consequence_counts = df["consequence"].value_counts().to_dict()

    def safe_extract(column):
        return df[column].dropna().tolist() if column in df.columns else []

    scores = {
        "sift4g_score": safe_extract("sift4g_score"),
        "polyphen2_hdiv_score": safe_extract("polyphen2_hdiv_score"),
        "cadd_phred": safe_extract("cadd_phred"),
    }

    return {
        "n_mutations": len(df),
        "consequences": consequence_counts,
        "score_counts": {k: len(v) for k, v in scores.items()},
        "score_means": {k: round(np.mean(v), 3) if v else None for k, v in scores.items()},
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--signature", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--rep", type=int, required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--outdir", default="results")
    parser.add_argument("--sigfile", default="signatures/COSMIC_v3.4_SBS_GRCh38.txt")
    parser.add_argument("--context", default="indices/context96.pkl")
    parser.add_argument("--exome", default="parquet")
    parser.add_argument("--fasta", default="ref/GRCh38.primary_assembly.genome.fa")
    parser.add_argument("--dbnsfp", default="ref/dbNSFP5.2a_grch38.gz")
    args = parser.parse_args()

    sig = args.signature
    outbase = Path(args.outdir) / sig / f"n_{args.n}"
    outbase.mkdir(parents=True, exist_ok=True)

    parquet_path = outbase / f"rep_{args.rep:02d}.annot.parquet"
    json_path = outbase / f"rep_{args.rep:02d}.summary.json"
    seed = args.seed if args.seed is not None else args.rep + args.n

    cmd = [
        "python", "scripts/simulate_mutations_1.py",
        "--signatures", args.sigfile,
        "--signature-name", sig,
        "--n", str(args.n),
        "--exome-dir", args.exome,
        "--context-index", args.context,
        "--fasta", args.fasta,
        "--out", str(parquet_path),
        "--seed", str(seed),
    ]

    subprocess.run(cmd, check=True)

    # Summarize
    df_summary = pq.read_table(parquet_path).to_pandas()
    summary = summarize_output(df_summary)
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)
    print(f"\n✅ Finished: {parquet_path} ({summary['n_mutations']} mutations)")

if __name__ == "__main__":
    main()
