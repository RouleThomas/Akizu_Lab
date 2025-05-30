#!/usr/bin/env python3
"""
Robust annotate_damage2.py that avoids UnicodeDecodeError:
- Uses subprocess to call tabix and decode output manually
- Ensures UTF-8 decoding with fallback on decode errors
"""
import sys, argparse, json, subprocess
from pathlib import Path
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq

def _first_numeric(field: str) -> float | None:
    for part in field.split(';'):
        try:
            return float(part)
        except ValueError:
            continue
    return None

class DBNSFP:
    IDX_SIFT4G_SCORE = 50 - 1
    IDX_SIFT4G_PRED  = 52 - 1
    IDX_PP2_HDIV_SCORE = 53 - 1
    IDX_PP2_HDIV_PRED  = 55 - 1
    IDX_CADD_PHRED     = 145 - 1

    def __init__(self, path: Path):
        self.path = path

    def query(self, chrom: str, pos: int, ref: str, alt: str) -> dict:
        chroms = [chrom.lstrip("chr"), f"chr{chrom.lstrip('chr')}"]
        for cand in chroms:
            try:
                result = subprocess.run(
                    ["tabix", str(self.path), f"{cand}:{pos}-{pos}"],
                    capture_output=True, check=True
                )
                lines = result.stdout.decode("utf-8", errors="replace").strip().split('\n')
            except subprocess.CalledProcessError:
                continue

            for line in lines:
                fields = line.strip().split('\t')
                if len(fields) <= max(
                    self.IDX_SIFT4G_SCORE,
                    self.IDX_SIFT4G_PRED,
                    self.IDX_PP2_HDIV_SCORE,
                    self.IDX_PP2_HDIV_PRED,
                    self.IDX_CADD_PHRED,
                    3
                ):
                    continue

                if fields[2] == ref and fields[3] == alt:
                    return {
                        "sift4g_score": _first_numeric(fields[self.IDX_SIFT4G_SCORE]),
                        "sift4g_pred":  fields[self.IDX_SIFT4G_PRED] or None,
                        "polyphen2_hdiv_score": _first_numeric(fields[self.IDX_PP2_HDIV_SCORE]),
                        "polyphen2_hdiv_pred":  fields[self.IDX_PP2_HDIV_PRED] or None,
                        "cadd_phred": _first_numeric(fields[self.IDX_CADD_PHRED]),
                    }
        return {
            "sift4g_score": None,
            "sift4g_pred": None,
            "polyphen2_hdiv_score": None,
            "polyphen2_hdiv_pred": None,
            "cadd_phred": None,
        }

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--dbnsfp", required=True, help="bgzip+tabix dbNSFP file")
    p.add_argument("in_parquet", help="Input Parquet (must have chr,pos,ref_base,alt_base)")
    p.add_argument("out_parquet", help="Output Parquet path")
    args = p.parse_args()

    df = pq.read_table(args.in_parquet).to_pandas()
    for col in ("chr", "pos", "ref_base", "alt_base"):
        if col not in df.columns:
            sys.exit(f"ERROR: missing column {col} in {args.in_parquet}")

    db = DBNSFP(Path(args.dbnsfp))
    annotations = {k: [] for k in [
        "sift4g_score", "sift4g_pred",
        "polyphen2_hdiv_score", "polyphen2_hdiv_pred",
        "cadd_phred"]}

    for chrom, pos, r_idx, a in zip(df.chr, df.pos, df.ref_base, df.alt_base):
        ref = "ACGT"[int(r_idx)]
        result = db.query(str(chrom), int(pos), ref, str(a))
        for k in annotations:
            annotations[k].append(result[k])

    for k in annotations:
        df[k] = annotations[k]

    tbl = pa.Table.from_pandas(df)
    pq.write_table(tbl, args.out_parquet, compression="zstd")

    n = len(df)
    print(json.dumps({
        "output": args.out_parquet,
        "rows": n,
        "NA_rates": {
            k: round(df[k].isna().mean() * 100, 2) for k in annotations
        }
    }, indent=2))

if __name__ == "__main__":
    main()

