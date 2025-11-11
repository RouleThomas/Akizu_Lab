#!/usr/bin/env python3
"""
Summarize simulation outputs into a single TSV focused on predictor category percentages,
using SCORE thresholds only (no *_pred columns).

Thresholds:
  SIFT4G:      <0.05 → Damaging; >=0.05 → Benign
  PolyPhen2:   <0.15 → Benign; 0.15–<0.85 → Possibly damaging; >=0.85 → Damaging
  CADD PHRED:  <2 → Benign; 2–<10 → Uncertain; 10–<20 → Likely damaging; >=20 → Damaging
"""

import argparse, json
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq


def _percent(num, den):
    return None if den in (None, 0) or num is None else 100.0 * float(num) / float(den)


def _col(df, names):
    for n in names:
        if n in df.columns:
            return n
    return None


def _sift4g_percentages(df):
    score_col = _col(df, ["sift4g_score", "SIFT4G_score", "sift_score"])
    if score_col is None:
        return dict(sift4g_n=0, sift4g_damaging_pct=None, sift4g_benign_pct=None)
    sc = pd.to_numeric(df[score_col], errors="coerce").dropna()
    if sc.empty:
        return dict(sift4g_n=0, sift4g_damaging_pct=None, sift4g_benign_pct=None)
    den = len(sc)
    return dict(
        sift4g_n=den,
        sift4g_damaging_pct=_percent((sc < 0.05).sum(), den),
        sift4g_benign_pct=_percent((sc >= 0.05).sum(), den),
    )


def _polyphen2_percentages(df):
    score_col = _col(df, ["polyphen2_hdiv_score", "polyphen2_score", "PPH2_HDIV_score", "PPH2_score"])
    if score_col is None:
        return dict(
            polyphen2_n=0,
            polyphen2_damaging_pct=None,
            polyphen2_possibly_damaging_pct=None,
            polyphen2_benign_pct=None,
        )
    sc = pd.to_numeric(df[score_col], errors="coerce").dropna()
    if sc.empty:
        return dict(
            polyphen2_n=0,
            polyphen2_damaging_pct=None,
            polyphen2_possibly_damaging_pct=None,
            polyphen2_benign_pct=None,
        )
    den = len(sc)
    damaging = (sc >= 0.85).sum()
    possibly = ((sc >= 0.15) & (sc < 0.85)).sum()
    benign = (sc < 0.15).sum()
    return dict(
        polyphen2_n=den,
        polyphen2_damaging_pct=_percent(damaging, den),
        polyphen2_possibly_damaging_pct=_percent(possibly, den),
        polyphen2_benign_pct=_percent(benign, den),
    )


def _cadd_percentages(df):
    phred_col = _col(df, ["cadd_phred", "CADD_phred", "caddPHRED"])
    if phred_col is None:
        return dict(
            cadd_n=0,
            cadd_damaging_pct=None,
            cadd_likely_damaging_pct=None,
            cadd_uncertain_pct=None,
            cadd_benign_pct=None,
        )
    sc = pd.to_numeric(df[phred_col], errors="coerce").dropna()
    if sc.empty:
        return dict(
            cadd_n=0,
            cadd_damaging_pct=None,
            cadd_likely_damaging_pct=None,
            cadd_uncertain_pct=None,
            cadd_benign_pct=None,
        )
    den = len(sc)
    benign = (sc < 2).sum()
    uncertain = ((sc >= 2) & (sc < 10)).sum()
    likely = ((sc >= 10) & (sc < 20)).sum()
    damaging = (sc >= 20).sum()
    return dict(
        cadd_n=den,
        cadd_damaging_pct=_percent(damaging, den),
        cadd_likely_damaging_pct=_percent(likely, den),
        cadd_uncertain_pct=_percent(uncertain, den),
        cadd_benign_pct=_percent(benign, den),
    )


def _summarize_one_replicate(json_file: Path):
    try:
        rep = int(json_file.stem.replace("rep_", "").split(".")[0])
    except Exception:
        print(f"⚠️ Skipping malformed file: {json_file}")
        return None
    try:
        n = int(json_file.parent.name.split("_")[1])
    except Exception:
        print(f"⚠️ Skipping unexpected folder name: {json_file.parent}")
        return None

    n_total = None
    try:
        with open(json_file) as f:
            data = json.load(f)
        n_total = data.get("n_mutations") if isinstance(data.get("n_mutations"), int) else None
    except Exception:
        pass

    parquet_file = json_file.with_name(json_file.name.replace("summary.json", "annot.parquet"))
    try:
        df = pq.read_table(parquet_file).to_pandas()
    except Exception as e:
        print(f"⚠️ Missing/invalid parquet for {json_file}: {e}")
        return None

    out = dict(n=n, rep=rep, n_total=n_total)
    out.update(_sift4g_percentages(df))
    out.update(_polyphen2_percentages(df))
    out.update(_cadd_percentages(df))
    return out


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--output", required=True)
    args = ap.parse_args()

    rows = []
    for json_file in Path(args.results_dir).rglob("rep_*.summary.json"):
        row = _summarize_one_replicate(json_file)
        if row:
            rows.append(row)

    cols = [
        "n","rep","n_total",
        "sift4g_n","sift4g_damaging_pct","sift4g_benign_pct",
        "polyphen2_n","polyphen2_damaging_pct","polyphen2_possibly_damaging_pct","polyphen2_benign_pct",
        "cadd_n","cadd_damaging_pct","cadd_likely_damaging_pct","cadd_uncertain_pct","cadd_benign_pct",
    ]

    df_out = pd.DataFrame(rows, columns=cols)
    df_out.sort_values(["n", "rep"], inplace=True)
    df_out.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Saved: {args.output}")


if __name__ == "__main__":
    main()
