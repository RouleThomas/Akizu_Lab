#!/usr/bin/env python3
"""
Robust COSMIC-style DBS (doublet base substitution) context index builder:
- Converts all base positions to 78 DBS contexts (with strand normalization)
- Output: context78.pkl mapping 0–77 (and -1) → list of row indices
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
OUT = HOME / "indices" / "context78.pkl"
OUT.parent.mkdir(parents=True, exist_ok=True)

# --- COSMIC canonical 78 DBS context IDs ---
# From COSMIC v3.3 DBS78 definition
# --- COSMIC canonical 78 DBS context IDs (order from user list) ---
CONTEXTS_78 = [
    "AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
    "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
    "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
    "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
    "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
    "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
    "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
    "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
    "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
    "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG",
]
CONTEXT2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_78)}

BASES = np.array(list("ACGT"))
def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT", "TGCA"))[::-1]

# --- Load reference FASTA ---
fasta = pysam.FastaFile(str(FASTA))
context_index: dict[int, list[int]] = defaultdict(list)

# --- Main loop ---
row_counter = 0
for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    tbl = pq.read_table(parquet, columns=["chr", "pos", "ref_base"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()
    refs = BASES[tbl["ref_base"].to_numpy()]

    for chrom, pos, ref in zip(chrs, poss, refs):
        # fetch 2-bp window for doublet context
        doublet = fasta.fetch(chrom, pos - 1, pos + 1).upper()
        if len(doublet) != 2 or "N" in doublet:
            context_index[-1].append(row_counter)
            row_counter += 1
            continue

        matched = False
        # Canonical orientation: leftmost pyrimidine (C/T)
        if doublet[0] in "CT":
            ref_dinuc = doublet
            for alt1 in "ACGT":
                for alt2 in "ACGT":
                    if alt1 + alt2 != ref_dinuc:
                        ctx = f"{ref_dinuc}>{alt1+alt2}"
                        if ctx in CONTEXT2ID:
                            context_index[CONTEXT2ID[ctx]].append(row_counter)
                            matched = True
        else:
            # reverse complement normalization
            rc = revcomp(doublet)
            ref_dinuc = rc
            if rc[0] in "CT":
                for alt1 in "ACGT":
                    for alt2 in "ACGT":
                        alt_rc = revcomp(alt1 + alt2)
                        if alt_rc != ref_dinuc:
                            ctx = f"{ref_dinuc}>{alt_rc}"
                            if ctx in CONTEXT2ID:
                                context_index[CONTEXT2ID[ctx]].append(row_counter)
                                matched = True

        if not matched:
            context_index[-1].append(row_counter)
        row_counter += 1

    print(f"Processed {parquet.name}, cumulative rows: {row_counter:,}")

# --- Ensure all 78 contexts are present ---
for i in range(78):
    context_index[i] = context_index.get(i, [])

# --- Save ---
full_index = {i: np.array(rows, dtype=np.int32) for i, rows in context_index.items()}
with open(OUT, "wb") as f:
    pickle.dump(full_index, f, protocol=4)

print(f"\n✅ Saved {OUT}")
print(f"✔️  Contexts stored: {len(full_index)} (should be 79 including fallback)")
print(f"🧠 Total indices stored (redundant): {sum(len(v) for v in full_index.values()):,}")
print(f"🔍 Unique positions covered: {len(np.unique(np.concatenate(list(full_index.values())))):,}")
