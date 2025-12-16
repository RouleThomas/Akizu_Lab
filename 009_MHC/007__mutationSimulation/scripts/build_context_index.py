#!/usr/bin/env python3
"""
Robust COSMIC-style context index builder:
- Converts all base positions to 96 SBS contexts (with strand normalization)
- Stores all positions: if unmatched, adds to fallback context -1
- Output: context96.pkl mapping 0–95 (and -1) → list of row indices
"""

from __future__ import annotations
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pysam


# --- Setup paths ---
HOME = Path.cwd()
PARQ_DIR = HOME / "parquet"
FASTA = HOME / "ref" / "GRCh38.primary_assembly.genome.fa"
OUT = HOME / "indices" / "context96.pkl"
OUT.parent.mkdir(parents=True, exist_ok=True)

# --- COSMIC canonical 96 SBS context IDs ---
BASES = "ACGT"
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref in "CT"
    for alt in ("AGT" if ref == "C" else "ACG")
    for l in BASES
    for r in BASES
]
CONTEXT2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}

# --- Utility ---
INT2BASE = np.array(list("ACGT"))

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# --- Load reference FASTA ---
fasta = pysam.FastaFile(str(FASTA))
context_index: dict[int, list[int]] = defaultdict(list)

# --- Main loop: process each chr parquet ---
row_counter = 0
for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    tbl = pq.read_table(parquet, columns=["chr", "pos", "ref_base"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()
    refs = INT2BASE[tbl["ref_base"].to_numpy()]

    for chrom, pos, ref in zip(chrs, poss, refs):
        tri = fasta.fetch(chrom, pos - 2, pos + 1).upper()
        if len(tri) != 3 or "N" in tri:
            context_index[-1].append(row_counter)
            row_counter += 1
            continue

        if tri[1] != ref:
            context_index[-1].append(row_counter)
            row_counter += 1
            continue

        matched = False

        # Handle canonical strand (ref is C or T)
        if ref in "CT":
            l, r = tri[0], tri[2]
            for alt in ("AGT" if ref == "C" else "ACG"):
                ctx = f"{l}[{ref}>{alt}]{r}"
                ctx_id = CONTEXT2ID.get(ctx)
                if ctx_id is not None:
                    context_index[ctx_id].append(row_counter)
                    matched = True

        # Handle reverse strand (ref is A or G)
        else:
            tri_rc = revcomp(tri)
            ref_rc = tri_rc[1]
            if ref_rc in "CT":
                l, r = tri_rc[0], tri_rc[2]
                for alt in ("AGT" if ref_rc == "C" else "ACG"):
                    ctx = f"{l}[{ref_rc}>{alt}]{r}"
                    ctx_id = CONTEXT2ID.get(ctx)
                    if ctx_id is not None:
                        context_index[ctx_id].append(row_counter)
                        matched = True

        # Store unmatched positions in fallback
        if not matched:
            context_index[-1].append(row_counter)

        row_counter += 1

    print(f"Processed {parquet.name}, cumulative rows: {row_counter:,}")

# --- Ensure all 96 contexts are present ---
for i in range(96):
    context_index[i] = context_index.get(i, [])

# --- Save ---
full_index = {i: np.array(rows, dtype=np.int32) for i, rows in context_index.items()}

with open(OUT, "wb") as f:
    pickle.dump(full_index, f, protocol=4)

print(f"\n✅ Saved {OUT}")
print(f"✔️  Contexts stored: {len(full_index)} (should be 97 including fallback)")
print(f"🧠 Total indices stored (redundant): {sum(len(v) for v in full_index.values()):,}")
print(f"🔍 Unique positions covered: {len(np.unique(np.concatenate(list(full_index.values())))):,}")
