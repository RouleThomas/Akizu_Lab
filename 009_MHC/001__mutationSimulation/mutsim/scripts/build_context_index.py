# build_context_index.py  (fixed – guarantees all 96 IDs)
"""
Build indices/context96.pkl mapping each of the 96 canonical SBS contexts
(ID 0‑95) to the *row indices* (int32) in the concatenated human‑exome table.

Fix: ensure **all 96 IDs exist** even though many share the same 3‑base
trinucleotide pool.  We first collect rows by trinucleotide key (left+ref+right)
then replicate that pool for every SBS substitution class that uses that
context.
"""
from __future__ import annotations

import pickle
from collections import defaultdict
from pathlib import Path

import numpy as np
import pyarrow.parquet as pq
import pysam


HOME = Path("/scr1/users/roulet/Akizu_Lab/009_MHC/001__mutationSimulation")
PARQ_DIR = HOME / "mutsim" / "parquet"
FASTA = HOME / "mutsim" / "ref" / "GRCh38.primary_assembly.genome.fa"
OUT = HOME / "mutsim" / "indices" / "context96.pkl"
OUT.parent.mkdir(parents=True, exist_ok=True)

# Canonical 96‑context order (matches COSMIC and simulate_mutations.py)
CONTEXTS_96 = [f"{l}{ref}{r}" for ref in "CT" for alt in ("AGT" if ref == "C" else "ACG") for l in "ACGT" for r in "ACGT"]

INT2BASE = np.array(list("ACGT"))  # ref_base column encoding 0‑3

fasta = pysam.FastaFile(str(FASTA))
tri_rows: dict[str, list[int]] = defaultdict(list)  # trinucleotide → rows
row_counter = 0

for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    tbl = pq.read_table(parquet, columns=["chr", "pos", "ref_base"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()
    refs = INT2BASE[tbl["ref_base"].to_numpy()]

    for chrom, pos, ref in zip(chrs, poss, refs):
        tri = fasta.fetch(chrom, pos - 2, pos + 1).upper()  # len 3 (L M R)
        if len(tri) != 3 or "N" in tri or ref not in "CT":
            row_counter += 1
            continue
        if tri[1] != ref:  # internal sanity check – rare due to repetitive DNA
            row_counter += 1
            continue
        tri_rows[tri].append(row_counter)
        row_counter += 1
    print(f"processed {parquet.name}, cumulative rows={row_counter:,}")

# ------------------------------------------------------------------
# Expand to full 96‑ID dict
full_index: dict[int, np.ndarray] = {}
for ctx_id, tri in enumerate(CONTEXTS_96):
    full_index[ctx_id] = np.array(tri_rows[tri], dtype=np.int32)

with open(OUT, "wb") as fh:
    pickle.dump(full_index, fh, protocol=4)
print("Saved", OUT, "with", len(full_index), "contexts (all IDs present)")

