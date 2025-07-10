#!/usr/bin/env python3
"""
Summarize simulation outputs (JSON + optional .annot.parquet fallback) into a single TSV.
Adds damaging score fractions from Parquet:
- frac_sift4g_damaging
- frac_polyphen2_damaging
- frac_cadd_damaging
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
            "frac_stop": (df["consequence"] == "stop").mean(),
            "frac_sift4g_damaging": (df["sift4g_pred"] == "D").mean(),
            "frac_polyphen2_damaging": df["polyphen2_hdiv_pred"].isin(["D", "P"]).mean(),
            "frac_cadd_damaging": (df["cadd_phred"] >= 20).mean()
        }
    except Exception as e:
        print(f"⚠️ Failed to read Parquet: {parquet_path} — {e}")
        return {
            "n_total": None, "frac_syn": None, "frac_missense": None, "frac_stop": None,
            "frac_sift4g_damaging": None, "frac_polyphen2_damaging": None, "frac_cadd_damaging": None
        }

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-dir", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    print(f"🔍 Searching in: {args.results_dir}")
    json_files = list(Path(args.results_dir).rglob("rep_*.summary.json"))
    print(f"📄 Found {len(json_files)} summary files.")

    rows = []
    for json_file in json_files:
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
                "frac_stop": consequences.get("stop", 0) / data["n_mutations"],
            })

            # Damaging fractions must be computed from Parquet
            parquet_file = json_file.with_name(json_file.name.replace("summary.json", "annot.parquet"))
            row.update(fallback_from_parquet(parquet_file))

        for key in ["sift4g_score", "polyphen2_hdiv_score", "cadd_phred"]:
            row[f"{key}_mean"] = data["score_means"].get(key)
            row[f"{key}_n"] = data["score_counts"].get(key)

        row["signature"] = json_file.parent.parent.name
        print(f"✅ Parsed: {json_file} → signature: {row['signature']}, rep: {row['rep']}, n: {row['n']}")
        rows.append(row)

    df_out = pd.DataFrame(rows)
    df_out.sort_values(["signature", "n", "rep"], inplace=True)
    df_out.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Saved: {args.output}")

if __name__ == "__main__":
    main()
