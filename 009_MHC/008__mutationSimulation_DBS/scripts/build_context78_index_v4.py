#!/usr/bin/env python3
from __future__ import annotations

import argparse, pickle
from pathlib import Path
from collections import defaultdict

import numpy as np
import pandas as pd
import pysam

DNA_COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s: str) -> str:
    return s.translate(DNA_COMP)[::-1]

def main():
    ap = argparse.ArgumentParser(description="Build strand-aware DBS78 context index (coding-start anchored).")
    ap.add_argument("--exome", default="parquet", help="Directory with chr*.parquet (your CDS table shards)")
    ap.add_argument("--fasta", required=True, help="Reference FASTA (indexed .fai)")
    ap.add_argument("--contexts-list", required=True, help="Text file with DBS78 contexts (one per line, e.g. AC>CA)")
    ap.add_argument("--out", default="indices/context78.pkl", help="Output pickle path")
    ap.add_argument("--max-chroms", type=int, default=None, help="Debug: limit number of chrom shards")
    args = ap.parse_args()

    exome_dir = Path(args.exome)
    shards = sorted(exome_dir.glob("chr*.parquet"))
    if not shards:
        raise SystemExit(f"No chr*.parquet found in {exome_dir}")
    if args.max_chroms:
        shards = shards[:args.max_chroms]

    # load DBS78 list
    ctx_names = [ln.strip() for ln in open(args.contexts_list) if ln.strip()]
    # map REF 2mer -> context id(s). For indexing we only need REF side.
    ref2_to_ids = defaultdict(list)
    for i, ctx in enumerate(ctx_names):
        ref2, _ = ctx.split(">")
        ref2_to_ids[ref2].append(i)

    wanted_ref2 = set(ref2_to_ids.keys())

    fa = pysam.FastaFile(args.fasta)

    pools = defaultdict(list)  # context_id -> list[row_index_global]
    row_offset = 0

    for fp in shards:
        df = pd.read_parquet(fp, columns=["chr","pos","strand"])
        # ensure python types
        chroms = df["chr"].astype(str).to_numpy()
        pos1   = df["pos"].astype(int).to_numpy()     # 1-based
        strand = df["strand"].astype(int).to_numpy()  # +/-1

        for i in range(len(df)):
            chrom = chroms[i]
            p = int(pos1[i])
            st = int(strand[i])

            # compute genomic slice for the DBS start in CODING direction
            # + strand: (p, p+1) => fetch [p-1, p+1)
            # - strand: (p-1, p) genomic, then revcomp to coding => fetch [p-2, p)
            try:
                if st == 1:
                    ref2 = fa.fetch(chrom, p-1, p+1).upper()  # two bases
                else:
                    # need p-1 exists
                    if p <= 1:
                        continue
                    ref2g = fa.fetch(chrom, p-2, p).upper()   # (p-1, p) in genomic order
                    ref2 = revcomp(ref2g)
            except Exception:
                continue

            if len(ref2) != 2 or "N" in ref2:
                continue

            # Only index contexts that exist in DBS78 reference set
            if ref2 not in wanted_ref2:
                continue

            # Add this row index to *all* contexts that share this REF2
            # (in DBS78 each ref2 corresponds to multiple alt2 choices)
            global_row = row_offset + i
            for ctx_id in ref2_to_ids[ref2]:
                pools[ctx_id].append(global_row)

        row_offset += len(df)
        print(f"✅ indexed {fp.name} (rows so far: {row_offset:,})")

    # convert to numpy arrays
    out = {int(k): np.asarray(v, dtype=np.int64) for k,v in pools.items()}
    out[-1] = np.array([], dtype=np.int64)  # sentinel if you like

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with open(out_path, "wb") as f:
        pickle.dump(out, f, protocol=pickle.HIGHEST_PROTOCOL)

    # small summary
    nonempty = sum(len(v) > 0 for v in out.values() if isinstance(v, np.ndarray))
    print(f"\n✅ wrote {out_path}")
    print(f"Contexts with pools: {nonempty} / {len(ctx_names)}")
    # quick check: report size of first few pools
    for k in sorted([k for k in out.keys() if k != -1])[:5]:
        print("ctx_id", k, "pool_size", len(out[k]))

if __name__ == "__main__":
    main()
