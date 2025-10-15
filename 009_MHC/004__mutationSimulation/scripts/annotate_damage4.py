#!/usr/bin/env python3
"""
Annotate simulated mutations with scores and predictions from dbNSFP,
prioritizing VEP canonical transcripts. Also outputs the canonical transcript ID.

Outputs a compressed Parquet table with score, prediction, and transcript ID columns.
"""

import sys, argparse, json, subprocess
from pathlib import Path
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
from statistics import median, mode

# --- Helper functions ---

def _split_or_none(field: str) -> list[str] | None:
    return field.split(';') if field and field != '.' else None

def _get_canonical_index(canon_flags: list[str]) -> int | None:
    if canon_flags is None:
        return None
    try:
        return canon_flags.index('YES')
    except ValueError:
        return None

def _get_value_at_index(values: list[str] | None, idx: int) -> str | None:
    if values and 0 <= idx < len(values):
        val = values[idx]
        return val if val not in {'.', ''} else None
    return None

def _get_float(val: str | None) -> float | None:
    try:
        return float(val) if val is not None else None
    except ValueError:
        return None

def _median_float(values: list[str] | None) -> float | None:
    if not values:
        return None
    floats = [_get_float(v) for v in values if v not in {None, '.', ''}]
    return median(floats) if floats else None

def _mode_str(values: list[str] | None) -> str | None:
    if not values:
        return None
    cleaned = [v for v in values if v not in {'.', ''}]
    try:
        return mode(cleaned) if cleaned else None
    except:
        return cleaned[0] if cleaned else None

# --- DBNSFP Parser ---

class DBNSFP:
    IDX_TRANSCRIPTID     = 15 - 1
    IDX_VEP_CANONICAL    = 26 - 1
    IDX_SIFT4G_SCORE     = 50 - 1
    IDX_SIFT4G_PRED      = 52 - 1
    IDX_PP2_HDIV_SCORE   = 53 - 1
    IDX_PP2_HDIV_PRED    = 55 - 1
    IDX_CADD_PHRED       = 144 - 1 # here updated!

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
                if len(fields) < 4:
                    continue  # Not enough fields to check ref/alt

                if fields[2] == ref and fields[3] == alt:
                    if len(fields) <= max(
                        self.IDX_TRANSCRIPTID,
                        self.IDX_VEP_CANONICAL,
                        self.IDX_SIFT4G_SCORE,
                        self.IDX_SIFT4G_PRED,
                        self.IDX_PP2_HDIV_SCORE,
                        self.IDX_PP2_HDIV_PRED,
                        self.IDX_CADD_PHRED,
                    ):
                        continue

                    canon_flags    = _split_or_none(fields[self.IDX_VEP_CANONICAL])
                    transcript_ids = _split_or_none(fields[self.IDX_TRANSCRIPTID])
                    sift_scores    = _split_or_none(fields[self.IDX_SIFT4G_SCORE])
                    sift_preds     = _split_or_none(fields[self.IDX_SIFT4G_PRED])
                    pp2_scores     = _split_or_none(fields[self.IDX_PP2_HDIV_SCORE])
                    pp2_preds      = _split_or_none(fields[self.IDX_PP2_HDIV_PRED])

                    idx = _get_canonical_index(canon_flags)
                    canonical_id = _get_value_at_index(transcript_ids, idx) if idx is not None else None

                    return {
                        "sift4g_score": _get_float(_get_value_at_index(sift_scores, idx)) if idx is not None else _median_float(sift_scores),
                        "sift4g_pred": _get_value_at_index(sift_preds, idx) if idx is not None else _mode_str(sift_preds),
                        "polyphen2_hdiv_score": _get_float(_get_value_at_index(pp2_scores, idx)) if idx is not None else _median_float(pp2_scores),
                        "polyphen2_hdiv_pred": _get_value_at_index(pp2_preds, idx) if idx is not None else _mode_str(pp2_preds),
                        "cadd_phred": _get_float(fields[self.IDX_CADD_PHRED]),
                        "canonical_transcript_id": canonical_id
                    }

        return {
            "sift4g_score": None,
            "sift4g_pred": None,
            "polyphen2_hdiv_score": None,
            "polyphen2_hdiv_pred": None,
            "cadd_phred": None,
            "canonical_transcript_id": None
        }

# --- Main script ---

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
    annotations = {
        "polyphen2_hdiv_score": [],
        "polyphen2_hdiv_pred": [],
        "sift4g_score": [],
        "sift4g_pred": [],
        "cadd_phred": [],
        "canonical_transcript_id": []
    }

    for chrom, pos, r_idx, a in zip(df.chr, df.pos, df.ref_base, df.alt_base):
        ref = "ACGT"[int(r_idx)]
        result = db.query(str(chrom), int(pos), ref, str(a))
        for k in annotations:
            annotations[k].append(result[k])

    for k, values in annotations.items():
        df[k] = values

    pq.write_table(pa.Table.from_pandas(df), args.out_parquet, compression="zstd")

    print(json.dumps({
        "output": args.out_parquet,
        "rows": len(df),
        "NA_rates": {k: round(df[k].isna().mean() * 100, 2) for k in annotations},
        "missing_canonical_ids": int(df["canonical_transcript_id"].isna().sum())
    }, indent=2))

if __name__ == "__main__":
    main()

