#!/usr/bin/env python3
"""
Wrapper script for a single SBS signature simulation + annotation run.
Outputs:
- Parquet with annotations
- JSON with summary stats (uses STRICT stop-gained)
"""

import argparse, json, time, os, sys, re
from pathlib import Path
import pyarrow as pa, pyarrow.parquet as pq
import numpy as np
from simulate_mutations import ExomeSampler, load_signatures, get_signature_probs
from annotate_damage4 import DBNSFP

# ---------- helpers ---------------------------------------------------------

# recognize common effect strings (VEP/snpEff/ANNOVAR)
RE_STOP_GAINED = re.compile(r'\b(?:stop[_ ]?gained|stopgain|nonsense(?:[_ ]?mutation)?)\b', re.I)
RE_STOP_OTHER  = re.compile(r'\b(?:stop[_ ]?retained|stop[_ ]?lost|stoploss)\b', re.I)
RE_MISSENSE    = re.compile(r'\b(?:missense(?:[_ ]?variant)?|nonsynonymous(?:[_ ]?snv)?)\b', re.I)
RE_SYNONYMOUS  = re.compile(r'\b(?:synonymous(?:[_ ]?variant)?|synonymous(?:[_ ]?snv)?)\b', re.I)
RE_HGVSP_STAR  = re.compile(r'p\.[A-Za-z]{3}\d+\*|[*]$', re.I)  # crude: “p.Gln123*” or trailing *

# context feasibility guard: can this SBS96 *ever* yield a stop?
CTX_RE = re.compile(r'^([ACGT])\[([ACGT])>([ACGT])\]([ACGT])$')
def context_can_stop(sig: str) -> bool:
    m = CTX_RE.match(sig)
    if not m:
        return True
    L, REF, ALT, R = m.groups()
    pos1 = (ALT == 'T' and R in {'A','G'})
    pos2 = (L == 'T' and ALT in {'A','G'} and R in {'A','G'})
    pos3 = (ALT in {'A','G'} and L in {'A','G'})
    return pos1 or pos2 or pos3

def _series_has(colname, df):  # convenience
    return colname in df.columns

def strict_stop_mask(df):
    """
    Return boolean mask for TRUE stop-gained:
      - prefer explicit effect labels (stop_gained/stopgain/nonsense_mutation),
      - else derive from AA fields if present (aa_ref != '*' & aa_alt == '*'),
      - else try HGVSp/AAChange strings,
      - else fall back to exact 'consequence' == 'stop_gained' only.
    """
    # 1) explicit effect labels in any column
    cand_cols = [c for c in df.columns if re.search(r'(consequence|csq|ann|annotation|exonicfunc|aachange|effect)', c, re.I)]
    if cand_cols:
        eff = df[cand_cols].astype(str).agg(" ".join, axis=1)
        m = eff.str.contains(RE_STOP_GAINED, na=False) & ~eff.str.contains(RE_STOP_OTHER, na=False)
        if m.any():
            return m

    # 2) explicit AA columns
    aa_cols = {c.lower(): c for c in df.columns}
    if 'aa_ref' in aa_cols and 'aa_alt' in aa_cols:
        m = (df[aa_cols['aa_ref']].astype(str) != '*') & (df[aa_cols['aa_alt']].astype(str) == '*')
        if m.any():
            return m

    # 3) HGVSp-like strings
    for c in ['HGVSp', 'hgvsp', 'AAChange.refGene', 'AAChange', 'amino_acids', 'Protein_position']:
        if c in df.columns:
            s = df[c].astype(str)
            m = s.str.contains(RE_HGVSP_STAR, na=False)
            if m.any():
                return m


    # 4) last resort: a simple 'consequence' column with coarse labels
    if 'consequence' in df.columns:
        s = df['consequence'].astype(str).str.lower()
        # If the schema is coarse (synonymous/missense/stop), treat 'stop' as gained.
        return s.isin({'stop', 'stop_gained', 'stopgain', 'nonsense_mutation'})


def summarize_output(df, signature):
    """
    Build strict summary dict from per-variant df.
    """
    n_total = len(df)

    # strict masks
    stop_gained = strict_stop_mask(df)

    # basic missense / synonymous (best-effort; independent of stop mask)
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

    # clamp impossible contexts
    if not context_can_stop(signature):
        stop_gained[:] = False

    # counts
    n_syn = int(synonymous.sum())
    n_mis = int(missense.sum())
    n_stop = int(stop_gained.sum())

    # fractions (denominator = total variants; change if you want coding-only)
    denom = max(int(n_total), 1)
    score_cols = ["sift4g_score", "polyphen2_hdiv_score", "cadd_phred"]
    scores = {k: df[k].dropna().tolist() for k in score_cols if k in df.columns}

    return {
        "n_mutations": int(n_total),
        "consequences": {
            "synonymous": n_syn,
            "missense": n_mis,
            "stop_gained": n_stop,     # <-- strict
        },
        "score_counts": {k: len(v) for k, v in scores.items()},
        "score_means": {k: (round(float(np.mean(v)), 3) if v else None) for k, v in scores.items()},
    }

# ---------- main ------------------------------------------------------------

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
    parser.add_argument("--dbnsfp", default="ref/dbNSFP5.2a_grch38.gz")
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

    sig_df = load_signatures(args.sigfile)
    probs = get_signature_probs(sig_df, sig)

    sampler = ExomeSampler(args.exome, args.context)
    if -1 in sampler.context_index:
        del sampler.context_index[-1]

    sampled_tbl, ctx_ids = sampler.sample(probs, args.n, rng)
    df = sampler.annotate(sampled_tbl, ctx_ids)

    df["pos_0based"] = df["pos"].astype(int)
    df["pos"] = df["pos_0based"] + 1   # now parquet pos is 1-based (VCF-style)

    # annotate dbNSFP scores (unchanged)
    db = DBNSFP(Path(args.dbnsfp))
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

    pq.write_table(pa.Table.from_pandas(df), parquet_path, compression="zstd")

    # ---------- STRICT summary here ----------
    summary = summarize_output(df, sig)
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"✅ Finished: {parquet_path} ({len(df)} mutations)")

if __name__ == "__main__":
    main()
