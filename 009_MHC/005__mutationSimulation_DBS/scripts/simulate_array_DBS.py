#!/usr/bin/env python3
"""
DBS-only simulator (no annotation).
Writes:
  <outdir>/<signature>/n_<n>/rep_<rep>.sim.parquet
  <outdir>/<signature>/n_<n>/rep_<rep>.sim.summary.json
  <outdir>/<signature>/n_<n>/rep_<rep>.index_debug.json
"""

import argparse, json, sys
from pathlib import Path
import numpy as np
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq

from simulate_mutations_DBS import ExomeSampler, load_signatures, get_signature_probs

def ensure_exists(label, path):
    if not Path(path).exists():
        sys.exit(f"❌ ERROR: {label} not found at: {path}")

def main():
    p = argparse.ArgumentParser()
    p.add_argument("--signature", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--rep", type=int, required=True)
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--outdir", default="results_contexts_DBS")
    p.add_argument("--sigfile", default="signatures/context_sigs_fixed_DBS.txt")
    p.add_argument("--context", default="indices/context78.pkl")
    p.add_argument("--exome", default="parquet")
    p.add_argument("--contexts_list", default="signatures/context_signature_list_DBS.txt")
    args = p.parse_args()

    # sanity
    ensure_exists("signature file", args.sigfile)
    ensure_exists("DBS context index", args.context)
    ensure_exists("exome dir", args.exome)
    ensure_exists("contexts list", args.contexts_list)

    sig = args.signature
    outbase = Path(args.outdir) / sig / f"n_{args.n}"
    outbase.mkdir(parents=True, exist_ok=True)
    parquet_path = outbase / f"rep_{args.rep:02d}.sim.parquet"
    json_path    = outbase / f"rep_{args.rep:02d}.sim.summary.json"
    debug_path   = outbase / f"rep_{args.rep:02d}.index_debug.json"

    seed = args.seed if args.seed is not None else args.rep + args.n
    rng = np.random.default_rng(seed)

    # context names order (used to map string-keyed indices)
    with open(args.contexts_list, "r") as fh:
        context_names = [ln.strip() for ln in fh if ln.strip()]

    # signatures aligned to that order
    sig_df = load_signatures(args.sigfile, context_order_path=args.contexts_list)
    probs  = get_signature_probs(sig_df, args.signature)

    sampler = ExomeSampler(args.exome, args.context, context_names=context_names)
    sampler.dump_debug(debug_path, probs_len=len(probs))

    # sample
    sampled_df, ctx_ids = sampler.sample(probs, args.n, rng)
    sampled_df["context_id"] = np.asarray(ctx_ids, dtype=int)  # positions into probs/order list
    sampled_df["signature"]  = sig
    sampled_df["rep"]        = int(args.rep)
    sampled_df["seed"]       = int(seed)
    sampled_df["n_requested"]= int(args.n)
    sampled_df["n_returned"] = len(sampled_df)

    pq.write_table(pa.Table.from_pandas(sampled_df), parquet_path, compression="zstd")

    ctx_counts = (
        pd.Series(sampled_df["context_id"]).value_counts().sort_index()
        if len(sampled_df) else pd.Series([], dtype=int)
    )
    if "chr" in sampled_df.columns:
        chr_counts = pd.Series(sampled_df["chr"]).value_counts().to_dict()
    elif "chrom" in sampled_df.columns:
        chr_counts = pd.Series(sampled_df["chrom"]).value_counts().to_dict()
    else:
        chr_counts = None

    summary = {
        "signature": sig,
        "rep": int(args.rep),
        "seed": int(seed),
        "n_requested": int(args.n),
        "n_returned": int(len(sampled_df)),
        "counts_by_context_id": {int(k): int(v) for k, v in ctx_counts.items()},
        "counts_by_chromosome": chr_counts,
        "paths": {
            "parquet": str(parquet_path),
            "index_debug": str(debug_path)
        }
    }
    with open(json_path, "w") as f:
        json.dump(summary, f, indent=2)

    print(f"✅ DBS simulation done: {parquet_path}  (n={len(sampled_df)})")

if __name__ == "__main__":
    main()
