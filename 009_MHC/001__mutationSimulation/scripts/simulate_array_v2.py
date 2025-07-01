#!/usr/bin/env python3
"""
Wrapper script for a single SBS signature simulation + annotation run.
Takes 3 arguments: signature name, number of mutations, and replicate ID.
Outputs:
- Annotated Parquet with SIFT/PolyPhen/CADD
- JSON with summary stats: mutation types, score distribution
"""

import argparse, json, time, os, sys
from pathlib import Path
import pyarrow as pa, pyarrow.parquet as pq
import numpy as np
import subprocess
from simulate_mutations import ExomeSampler, load_signatures, get_signature_probs
from annotate_damage2 import DBNSFP

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
    parser.add_argument("--signature", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--rep", type=int, required=True)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--outdir", default="/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/results") # MODIFIED
    parser.add_argument("--sigfile", default="/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/signatures/COSMIC_v3.4_SBS_GRCh38.txt") # MODIFIED
    parser.add_argument("--context", default="/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/indices/context96.pkl") # MODIFIED
    parser.add_argument("--exome", default="/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/parquet") # MODIFIED
    parser.add_argument("--dbnsfp", default="/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation/mutsim/ref/dbNSFP5.1a_grch38.gz") # MODIFIED
    args = parser.parse_args()

    # Fail-fast path check
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
    sampled_tbl, ctx_ids = sampler.sample(probs, args.n, rng)
    df = sampler.annotate(sampled_tbl, ctx_ids)

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
    with open(json_path, "w") as f:
        json.dump(summarize_output(df), f, indent=2)

    print(f"✅ Finished: {parquet_path} ({len(df)} mutations)")

if __name__ == "__main__":
    main()

