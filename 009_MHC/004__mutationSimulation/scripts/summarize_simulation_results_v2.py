#!/usr/bin/env python3
"""
Summarize simulation outputs into a single TSV with strict stop-gained.

- Prefer variant-level parquet (rep_*.annot.parquet) **only if** coding/effect
  annotations are detected and parsed.
- Otherwise, compute strict stop-gained from JSON:
  use 'stop_gained' if present, else 'stop - stop_retained - stop_lost'.
- Add --source {auto,json,parquet} to force behavior, and --debug for logging.
"""

import argparse, json, re, os
from pathlib import Path
import pandas as pd
import pyarrow.parquet as pq

# ---------------- regexes (case-insensitive) ----------------
RE_STOP_GAINED = re.compile(r'\b(?:stop[_ ]?gained|stopgain|nonsense(?:[_ ]?mutation)?)\b', re.I)
RE_STOP_OTHER  = re.compile(r'\b(?:stop[_ ]?retained|stop[_ ]?lost|stoploss)\b', re.I)
RE_MISSENSE    = re.compile(r'\b(?:missense(?:[_ ]?variant)?|nonsynonymous(?:[_ ]?snv)?)\b', re.I)
RE_SYNONYMOUS  = re.compile(r'\b(?:synonymous(?:[_ ]?variant)?|synonymous(?:[_ ]?snv)?)\b', re.I)
RE_CODING_ANY  = re.compile(
    r'\b(?:synonymous|missense|nonsynonymous|stop[_ ]?gained|stopgain|nonsense|'
    r'stop[_ ]?lost|stoploss|stop[_ ]?retained|start[_ ]?lost|frameshift|inframe|coding)\b', re.I
)

VEPCOLS = re.compile(r'(consequence|csq|ann|annotation|exonicfunc|aachange|hgvsp|amino|info)', re.I)

def _log(debug: bool, *msg):
    if debug:
        print(*msg)

def _pick_effect_cols(df: pd.DataFrame):
    return [c for c in df.columns if VEPCOLS.search(c)]

def _collapse_to_variants(df: pd.DataFrame) -> pd.DataFrame:
    # try to deduplicate transcript-level rows to variant-level rows
    keys = []
    for k in df.columns:
        lk = k.lower()
        if lk in {"chrom","chr","contig","pos","position","start","ref","alt","ref_allele","alt_allele","end"}:
            keys.append(k)
    if {"ref","alt"} <= {k.lower() for k in keys} and any(k.lower() in {"pos","position","start"} for k in keys):
        return df.drop_duplicates(subset=keys)
    return df

def _counts_from_parquet(parquet_path: Path, debug: bool):
    """Return dict or None if parquet unusable."""
    try:
        df = pq.read_table(parquet_path).to_pandas()
    except Exception:
        _log(debug, f"  parquet not found/failed: {parquet_path}")
        return None
    if df.empty:
        _log(debug, f"  parquet empty: {parquet_path}")
        return None

    df = _collapse_to_variants(df)
    # Case 1: we have a single 'consequence' column with simple labels
    for cname in df.columns:
        if cname.lower() == "consequence":
            s = df[cname].astype(str).str.lower()
            coding_mask = s.str.contains(RE_CODING_ANY, na=False) | s.isin({"synonymous","missense","stop","stop_gained","stopgain","nonsense_mutation"})
            if coding_mask.sum() == 0:
                break
            stop_gained = s.isin({"stop_gained","stopgain","nonsense_mutation","stop"})  # 'stop' treated as gained here
            synonymous  = s.isin({"synonymous","synonymous_variant"})
            missense    = s.isin({"missense","missense_variant","nonsynonymous","nonsynonymous snv"})
            denom = max(int(coding_mask.sum()), 1)
            _log(debug, f"  parquet used (simple consequence col='{cname}') coding={coding_mask.sum()} stop={int(stop_gained.sum())}")
            return {
                "n_total": int(coding_mask.sum()),
                "frac_syn": float((synonymous & coding_mask).sum() / denom),
                "frac_missense": float((missense & coding_mask).sum() / denom),
                "frac_stop": float((stop_gained & coding_mask).sum() / denom),
            }

    # Case 2: glue likely effect columns and regex them
    effect_cols = _pick_effect_cols(df)
    if not effect_cols:
        _log(debug, "  parquet has no effect-like columns; skipping")
        return None

    eff = df[effect_cols].astype(str).agg(" ".join, axis=1).str.lower()
    coding_mask = eff.str.contains(RE_CODING_ANY, na=False)
    if coding_mask.sum() == 0:
        _log(debug, "  parquet effect strings found, but no coding keywords; skipping")
        return None

    stop_gained = eff.str.contains(RE_STOP_GAINED, na=False) & ~eff.str.contains(RE_STOP_OTHER, na=False)
    synonymous  = eff.str.contains(RE_SYNONYMOUS, na=False) & ~stop_gained
    missense    = eff.str.contains(RE_MISSENSE,   na=False) & ~stop_gained
    denom = max(int(coding_mask.sum()), 1)
    _log(debug, f"  parquet used (regex) coding={coding_mask.sum()} stop={int((stop_gained & coding_mask).sum())}")
    return {
        "n_total": int(coding_mask.sum()),
        "frac_syn": float((synonymous & coding_mask).sum() / denom),
        "frac_missense": float((missense   & coding_mask).sum() / denom),
        "frac_stop": float((stop_gained    & coding_mask).sum() / denom),
    }

def _norm_keys(d: dict) -> dict:
    out = {}
    for k, v in d.items():
        k2 = k.lower().replace(" ", "_").replace("-", "_")
        out[k2] = v
    return out

def _counts_from_json(data: dict, debug: bool):
    cons = _norm_keys(data.get("consequences", {}))
    nmut = max(int(data.get("n_mutations", 0)), 1)
    if "stop_gained" in cons:
        stop_strict = cons["stop_gained"]
    else:
        stop_strict = cons.get("stop", 0) - cons.get("stop_retained", 0) - cons.get("stop_lost", 0) - cons.get("stoploss", 0)
        stop_strict = max(stop_strict, 0)
    _log(debug, f"  json used coding≈{nmut} stop_strict={stop_strict}")
    return {
        "n_total": int(data.get("n_mutations", 0)),
        "frac_syn": cons.get("synonymous", 0) / nmut,
        "frac_missense": cons.get("missense", 0) / nmut,
        "frac_stop": stop_strict / nmut,
    }

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--results-dir", required=True)
    ap.add_argument("--output", required=True)
    ap.add_argument("--source", choices=["auto","json","parquet"], default="auto",
                    help="Where to take effects from (default: auto)")
    ap.add_argument("--debug", action="store_true", help="Verbose logging per replicate")
    args = ap.parse_args()

    rows = []
    for json_file in Path(args.results_dir).rglob("rep_*.summary.json"):
        try:
            rep = int(json_file.stem.replace("rep_", "").split(".")[0])
        except ValueError:
            print(f"⚠️  Skipping malformed file: {json_file}")
            continue

        parent = json_file.parent
        try:
            n = int(parent.name.split("_")[1])
        except Exception:
            print(f"⚠️  Skipping folder with unexpected name: {parent}")
            continue

        with open(json_file) as f:
            data = json.load(f)

        parquet_file = json_file.with_name(json_file.name.replace("summary.json", "annot.parquet"))
        if args.debug:
            print(f"• {json_file}")

        if args.source == "parquet":
            counts = _counts_from_parquet(parquet_file, args.debug)
            if counts is None:
                raise SystemExit("Requested --source parquet but no usable annotations found.")
        elif args.source == "json":
            counts = _counts_from_json(data, args.debug)
        else:
            # auto: try parquet; if unusable, use JSON
            counts = _counts_from_parquet(parquet_file, args.debug)
            if counts is None:
                counts = _counts_from_json(data, args.debug)

        score_means = data.get("score_means", {})
        score_counts = data.get("score_counts", {})
        rows.append({
            "n": n,
            "rep": rep,
            "n_total": counts["n_total"],
            "frac_syn": counts["frac_syn"],
            "frac_missense": counts["frac_missense"],
            "frac_stop": counts["frac_stop"],   # strict stop-gained
            "sift4g_score_mean": score_means.get("sift4g_score"),
            "sift4g_score_n": score_counts.get("sift4g_score"),
            "polyphen2_hdiv_score_mean": score_means.get("polyphen2_hdiv_score"),
            "polyphen2_hdiv_score_n": score_counts.get("polyphen2_hdiv_score"),
            "cadd_phred_mean": score_means.get("cadd_phred"),
            "cadd_phred_n": score_counts.get("cadd_phred"),
        })

    df_out = pd.DataFrame(rows).sort_values(["n", "rep"])
    df_out.to_csv(args.output, sep="\t", index=False)
    print(f"✅ Saved: {args.output}")

if __name__ == "__main__":
    main()
