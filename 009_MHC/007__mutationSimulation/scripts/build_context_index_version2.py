#!/usr/bin/env python3
"""
COSMIC-style 96 SBS context index builder (strand-normalized).

Input:
  parquet/chr*.parquet with columns: chr, pos (1-based), ref_base (0:A 1:C 2:G 3:T 4:N)

Output:
  indices/context96.pkl mapping:
    0..95  -> np.int32 array of global row indices
    -1     -> fallback (unusable / mismatch / N / boundary) global row indices
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
INT2BASE = np.array(list("ACGTN"))

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# --- Load reference FASTA ---
fasta = pysam.FastaFile(str(FASTA))
context_index: dict[int, list[int]] = defaultdict(list)

# Cache chromosome lengths
chr_len_cache: dict[str, int] = {}

row_counter = 0

for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    tbl = pq.read_table(parquet, columns=["chr", "pos", "ref_base"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()
    refs = INT2BASE[tbl["ref_base"].to_numpy()]

    for chrom, pos, ref in zip(chrs, poss, refs):
        if chrom not in chr_len_cache:
            chr_len_cache[chrom] = fasta.get_reference_length(chrom)
        chrom_len = chr_len_cache[chrom]

        # Boundary guard: need pos-1 and pos+1 to exist (1-based)
        if pos < 2 or pos > chrom_len - 1:
            context_index[-1].append(row_counter)
            row_counter += 1
            continue

        # pysam.fetch uses 0-based, end-exclusive
        # Want tri centered on pos (1-based): [pos-1, pos, pos+1] => start0=pos-2, end0=pos+1
        tri = fasta.fetch(chrom, int(pos) - 2, int(pos) + 1).upper()

        if len(tri) != 3 or "N" in tri:
            context_index[-1].append(row_counter)
            row_counter += 1
            continue

        if ref == "N" or tri[1] != ref:
            context_index[-1].append(row_counter)
            row_counter += 1
            continue

        matched = False

        # Canonical strand (ref is C or T)
        if ref in "CT":
            l, r = tri[0], tri[2]
            for alt in ("AGT" if ref == "C" else "ACG"):
                ctx_id = CONTEXT2ID.get(f"{l}[{ref}>{alt}]{r}")
                if ctx_id is not None:
                    context_index[ctx_id].append(row_counter)
                    matched = True
        else:
            # Reverse strand: normalize to C/T by reverse-complementing tri
            tri_rc = revcomp(tri)
            ref_rc = tri_rc[1]
            if ref_rc in "CT":
                l, r = tri_rc[0], tri_rc[2]
                for alt in ("AGT" if ref_rc == "C" else "ACG"):
                    ctx_id = CONTEXT2ID.get(f"{l}[{ref_rc}>{alt}]{r}")
                    if ctx_id is not None:
                        context_index[ctx_id].append(row_counter)
                        matched = True

        if not matched:
            context_index[-1].append(row_counter)

        row_counter += 1

    print(f"Processed {parquet.name}, cumulative rows: {row_counter:,}")

# Ensure all 96 contexts exist
for i in range(96):
    context_index[i] = context_index.get(i, [])

# Save as numpy arrays
full_index = {i: np.array(rows, dtype=np.int32) for i, rows in context_index.items()}

with open(OUT, "wb") as f:
    pickle.dump(full_index, f, protocol=4)

# Robust summary printing
all_arrays = [v for v in full_index.values() if len(v) > 0]
unique_n = int(len(np.unique(np.concatenate(all_arrays)))) if all_arrays else 0

print(f"\n✅ Saved {OUT}")
print(f"✔️  Contexts stored: {len(full_index)} (should be 97 including fallback)")
print(f"🧠 Total indices stored (redundant): {sum(len(v) for v in full_index.values()):,}")
print(f"🔍 Unique positions covered: {unique_n:,}")
