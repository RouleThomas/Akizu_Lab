#!/usr/bin/env python3
"""
Wrapper script for a single SBS signature simulation + dbNSFP annotation run.
Assumes:
  - exome parquet files are 1-based positions (pos matches IGV/VCF style)
  - context96.pkl was built using those same 1-based positions
Outputs:
  - Parquet with dbNSFP annotations
  - JSON with summary stats (STRICT stop-gained/stop-lost best-effort)
"""

import argparse, json, sys, re
from pathlib import Path

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq

from simulate_mutations_v2 import ExomeSampler, load_signatures, get_signature_probs
from annotate_damage5 import DBNSFP


# ---------- helpers ---------------------------------------------------------

RE_STOP_GAINED = re.compile(r'\b(?:stop[_ ]?gained|stopgain|nonsense(?:[_ ]?mutation)?)\b', re.I)
RE_STOP_OTHER  = re.compile(r'\b(?:stop[_ ]?retained|stop[_ ]?lost|stoploss)\b', re.I)
RE_MISSENSE    = re.compile(r'\b(?:missense(?:[_ ]?variant)?|nonsynonymous(?:[_ ]?snv)?)\b', re.I)
RE_SYNONYMOUS  = re.compile(r'\b(?:synonymous(?:[_ ]?variant)?|synonymous(?:[_ ]?snv)?)\b', re.I)
RE_HGVSP_STAR  = re.compile(r'p\.[A-Za-z]{3}\d+\*|[*]$', re.I)

RE_STOP_LOST = re.compile(r'\b(?:stop[_ ]?lost|stoploss)\b', re.I)

def strict_stop_mask(df):
    # 1) explicit effect labels
    cand_cols = [c for c in df.columns if re.search(r'(consequence|csq|ann|annotation|exonicfunc|aachange|effect)', c, re.I)]
    if cand_cols:
        eff = df[cand_cols].astype(str).agg(" ".join, axis=1)
        m = eff.str.contains(RE_STOP_GAINED, na=False) & ~eff.str.contains(RE_STOP_OTHER, na=False)
        if m.any():
            return m

    # 2) AA columns if present
    aa_cols = {c.lower(): c for c in df.columns}
    if 'aa_ref' in aa_cols and 'aa_alt' in aa_cols:
        m = (df[aa_cols['aa_ref']].astype(str) != '*') & (df[aa_cols['aa_alt']].astype(str) == '*')
        if m.any():
            return m

    # 3) HGVSp-like strings
    for c in ['HGVSp', 'hgvsp', 'AAChange.refGene', 'AAChange', 'amino_acids']:
        if c in df.columns:
            s = df[c].astype(str)
            m = s.str.contains(RE_HGVSP_STAR, na=False)
            if m.any():
                return m

    # 4) coarse consequence fallback
    if 'consequence' in df.columns:
        s = df['consequence'].astype(str).str.lower()
        return s.isin({'stop', 'stop_gained', 'stopgain', 'nonsense_mutation'})

    return np.zeros(len(df), dtype=bool)

def strict_stop_lost_mask(df):
    cand_cols = [c for c in df.columns if re.search(r'(consequence|csq|ann|annotation|exonicfunc|aachange|effect)', c, re.I)]
    if cand_cols:
        eff = df[cand_cols].astype(str).agg(" ".join, axis=1)
        m = eff.str.contains(RE_STOP_LOST, na=False) & ~eff.str.contains(RE_STOP_GAINED, na=False)
        if m.any():
            return m

    aa_cols = {c.lower(): c for c in df.columns}
    if 'aa_ref' in aa_cols and 'aa_alt' in aa_cols:
        m = (df[aa_cols['aa_ref']].astype(str) == '*') & (df[aa_cols['aa_alt']].astype(str) != '*')
        if m.any():
            return m

    if 'consequence' in df.columns:
        s = df['consequence'].astype(str).str.lower()
        return s.isin({'stop_lost', 'stoploss', 'stop_lost_variant'})

    return np.zeros(len(df), dtype=bool)

def summarize_output(df):
    n_total = int(len(df))
    stop_gained = strict_stop_mask(df)
    stop_lost   = strict_stop_lost_mask(df)
    stop_lost   = stop_lost & ~stop_gained

    missense = np.zeros(n_total, dtype=bool)
    synonymous = np.zeros(n_total, dtype=bool)

    if 'consequence' in df.columns:
        s = df['consequence'].astype(str).str.lower()
        missense |= s.isin(['missense', 'missense_variant', 'nonsynonymous', 'nonsynonymous snv'])
        synonymous |= s.isin(['synonymous', 'synonymous_variant', 'synonymous snv'])
    else:
        cand_cols = [c for c in df.columns if re.search(r'(consequence|csq|ann|annotation|exonicfunc|aachange|effect)', c, re.I)]
        if cand_cols:
            eff = df[cand_cols].astype(str).agg(" ".join, axis=1)
            missense |= eff.str.contains(RE_MISSENSE, na=False)
            synonymous |= eff.str.contains(RE_SYNONYMOUS, na=False)

    score_cols = ["sift4g_score", "polyphen2_hdiv_score", "cadd_phred"]
    scores = {k: df[k].dropna().tolist() for k in score_cols if k in df.columns}

    return {
        "n_mutations": n_total,
        "consequences": {
            "synonymous": int(synonymous.sum()),
            "missense": int(missense.sum()),
            "stop_gained": int(stop_gained.sum()),
            "stop_lost": int(stop_lost.sum()),
        },
        "score_counts": {k: len(v) for k, v in scores.items()},
        "score_means": {k: (round(float(np.mean(v)), 3) if v else None) for k, v in scores.items()},
    }

def normalize_chrom_for_dbnsfp(chrom: str) -> str:
    # Try to be tolerant: DBNSFP sometimes uses "1" not "chr1"
    c = str(chrom)
    return c[3:] if c.lower().startswith("chr") else c


# ---------- main ------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(allow_abbrev=False)
    parser.add_argument("--signature", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--rep", type=int, required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--outdir", default="results_v5")
    parser.add_argument("--sigfile", default="signatures/COSMIC_with_flat.txt")
    parser.add_argument("--context", default="indices/context96.pkl")
    parser.add_argument("--exome", default="parquet")
    parser.add_argument("--dbnsfp", default="ref/dbNSFP5.3a_grch38.gz")
    args = parser.parse_args()

    for label, path in {
        "signature file": args.sigfile,
        "context index": args.context,
        "exome dir": args.exome,
        "dbNSFP file": args.dbnsfp
    }.items():
        if not Path(path).exists():
            sys.exit(f"❌ ERROR: {label} not found at: {path}")

    sig = args.signature
    outbase = Path(args.outdir) / sig / f"n_{args.n}"
    outbase.mkdir(parents=True, exist_ok=True)
    parquet_path = outbase / f"rep_{args.rep:02d}.annot.parquet"
    json_path = outbase / f"rep_{args.rep:02d}.summary.json"

    seed = args.seed if args.seed is not None else args.rep + args.n
    rng = np.random.default_rng(seed)

    # Load signature probabilities (SBS96 order must match your signature file format)
    sig_df = load_signatures(args.sigfile)
    probs = get_signature_probs(sig_df, sig)

    sampler = ExomeSampler(args.exome, args.context)
    # drop fallback bucket if present
    if hasattr(sampler, "context_index") and (-1 in sampler.context_index):
        del sampler.context_index[-1]

    sampled_tbl, ctx_ids = sampler.sample(probs, args.n, rng)
    df = sampler.annotate(sampled_tbl, ctx_ids)

    # IMPORTANT: with your new parquet, df["pos"] is already 1-based. DO NOT shift it.
    df["pos"] = df["pos"].astype(int)

    # dbNSFP annotate
    db = DBNSFP(Path(args.dbnsfp))
    annotations = {k: [] for k in [
        "polyphen2_hdiv_score", "polyphen2_hdiv_pred",
        "sift4g_score", "sift4g_pred", "cadd_phred"
    ]}

    for chrom, pos, r_idx, alt in zip(df["chr"], df["pos"], df["ref_base"], df["alt_base"]):
        # guard against ref_base == 4 (N)
        ri = int(r_idx)
        if ri < 0 or ri > 3:
            result = {k: None for k in annotations}
        else:
            ref = "ACGT"[ri]
            chrom_q = normalize_chrom_for_dbnsfp(chrom)
            # Try normalized chrom; if your DBNSFP class expects chr-prefix, adjust there instead.
            result = db.query(str(chrom_q), int(pos), ref, str(alt))

        for key in annotations:
            annotations[key].append(result.get(key))

    for k, v in annotations.items():
        df[k] = v

    pq.write_table(pa.Table.from_pandas(df), parquet_path, compression="zstd")

    summary = summarize_output(df)
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"✅ Finished: {parquet_path} ({len(df)} mutations)")

if __name__ == "__main__":
    main()
