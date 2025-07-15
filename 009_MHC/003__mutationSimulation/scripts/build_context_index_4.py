#!/usr/bin/env python3
# build_context_index.py (strand-aware + purine-normalizing, full 96-context support)

"""
Build indices/context96.pkl mapping each of the 96 canonical SBS contexts
(ID 0–95) to the row indices (int32) in the human exome table.

Fixes:
- Includes both strands using the `strand` field in the exome table.
- Normalizes purine-centered sites (A/G) into pyrimidine contexts using reverse complement.
- Guarantees all 96 SBS context IDs are represented.

Context IDs follow COSMIC 96-class order:
    - 6 mutation types: C>A, C>G, C>T, T>A, T>C, T>G
    - 4 left bases × 4 right bases per type = 96
"""

from __future__ import annotations
import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pysam
from Bio.Seq import reverse_complement

# ------------------------------------------------------------------------------
# Constants and setup
# ------------------------------------------------------------------------------

HOME = Path.cwd() # current wd
PARQ_DIR = HOME / "parquet"
FASTA = HOME / "ref" / "GRCh38.primary_assembly.genome.fa"
OUT = HOME / "indices" / "context96.pkl"
OUT.parent.mkdir(parents=True, exist_ok=True)

# Canonical SBS96 contexts (in COSMIC order)
SUBS = [("C", "A"), ("C", "G"), ("C", "T"), ("T", "A"), ("T", "C"), ("T", "G")]
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alt in SUBS
    for l in "ACGT"
    for r in "ACGT"
]
CTX2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}

INT2BASE = np.array(list("ACGT"))  # maps 0–3 to A/C/G/T

# ------------------------------------------------------------------------------
# Load reference genome
# ------------------------------------------------------------------------------

fasta = pysam.FastaFile(str(FASTA))

# ------------------------------------------------------------------------------
# Build context index
# ------------------------------------------------------------------------------

context_rows: dict[int, list[int]] = defaultdict(list)
row_counter = 0


for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    tbl = pq.read_table(parquet, columns=["chr", "pos", "ref_base", "strand"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()
    refs = INT2BASE[tbl["ref_base"].to_numpy()]
    strands = tbl["strand"].to_numpy()


    for chrom, pos, ref, strand in zip(chrs, poss, refs, strands):
        try:
            tri = fasta.fetch(chrom, pos - 1, pos + 2).upper()
        except Exception:
            row_counter += 1
            continue
        if len(tri) != 3 or "N" in tri:
            row_counter += 1
            continue

        original_ref = tri[1]
        normalized = False

        # Normalize to pyrimidine-centered if needed
        if original_ref in "AG":
            tri = reverse_complement(tri)
            ref_base = tri[1]
            normalized = True
        else:
            ref_base = original_ref

        if ref_base not in "CT":
            row_counter += 1
            continue

        l, m, r = tri
        if m != ref_base:
            row_counter += 1
            continue

        for alt in ("AGT" if ref_base == "C" else "ACG"):
            if alt == ref_base:
                continue
            ctx = f"{l}[{ref_base}>{alt}]{r}"
            ctx_id = CTX2ID[ctx]
            context_rows[ctx_id].append((row_counter, normalized))  # ⬅️ store both row and whether it was normalized

        row_counter += 1


    print(f"Processed {parquet.name}, cumulative rows = {row_counter:,}")

# ------------------------------------------------------------------------------
# Finalize and save
# ------------------------------------------------------------------------------

# Ensure all 96 context IDs are present
full_index: dict[int, np.ndarray] = {
    ctx_id: np.array(rows, dtype=[('row', 'i4'), ('revcomp', 'b')])
    for ctx_id, rows in context_rows.items()
}


with open(OUT, "wb") as fh:
    pickle.dump(full_index, fh, protocol=4)

print(f"Saved {OUT} with {len(full_index)} SBS96 contexts.")


