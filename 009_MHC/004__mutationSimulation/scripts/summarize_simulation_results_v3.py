#!/usr/bin/env python3
"""
Summarize simulation outputs (JSON + optional .annot.parquet fallback) into a single TSV.
"""

import argparse, json
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

def fallback_from_parquet(parquet_path):
    try:
        df = pq.read_table(parquet_path).to_pandas()
        return {
            "n_total": len(df),
            "frac_syn": (df["consequence"] == "synonymous").mean(),
            "frac_missense": (df["consequence"] == "missense").mean(),
            "frac_stop": (df["consequence"] == "stop_gained").mean()
        }
    except Exception:
        return {"n_total": None, "frac_syn": None, "frac_missense": None, "frac_stop": None}

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    rows = []
    for json_file in Path(args.results_dir).rglob("rep_*.summary.json"):
        try:
            rep_str = json_file.stem.replace("rep_", "").split(".")[0]
            rep = int(rep_str)
        except ValueError:
            print(f"⚠️ Skipping malformed file: {json_file}")
            continue

        parent = json_file.parent
        try:
            n = int(parent.name.split("_")[1])
        except Exception:
            print(f"⚠️ Skipping folder with unexpected name: {parent}")
            continue

        with open(json_file) as f:
            data = json.load(f)

        row = {"n": n, "rep": rep}
        fallback_needed = "consequences" not in data

        if fallback_needed:
            parquet_file = json_file.with_name(json_file.name.replace("summary.json", "annot.parquet"))
            row.update(fallback_from_parquet(parquet_file))
        else:
            consequences = data["consequences"]
            row.update({
                "n_total": data["n_mutations"],
                "frac_syn": consequences.get("synonymous", 0) / data["n_mutations"],
                "frac_missense": consequences.get("missense", 0) / data["n_mutations"],
                "frac_stop": consequences.get("stop_gained", 0) / data["n_mutations"],
            })

        for key in ["sift4g_score", "polyphen2_hdiv_score", "cadd_phred"]:
            row[f"{key}_mean"] = data["score_means"].get(key)
            row[f"{key}_n"] = data["score_counts"].get(key)

        rows.append(row)

    df_out = pd.DataFrame(rows)
    df_out.sort_values(["n", "rep"], inplace=True)
    df_out.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Saved: {args.output}")

if __name__ == "__main__":
    main()

