#!/usr/bin/env python3
"""
Build DBS context index (0..77 and -1) in the SAME representation
as your signature file (i.e., keep forward-strand ref doublets; NO
pyrimidine-left canonicalization).

Assumes:
  - parquet pos is 1-based
  - DBS context is anchored at the LEFT base of the dinucleotide at `pos`
"""

from __future__ import annotations
from pathlib import Path
from collections import defaultdict
import numpy as np
import pyarrow.parquet as pq
import pysam, pickle

HOME   = Path.cwd()
PARQ   = HOME / "parquet"   # chr*.parquet used by simulator
FASTA  = HOME / "ref" / "GRCh38.primary_assembly.genome.fa"
OUT    = HOME / "indices" / "context78.pkl"
OUT.parent.mkdir(parents=True, exist_ok=True)

# Use YOUR DBS78 list (same order and spelling as your signature file)
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
CTX2ID = {c: i for i, c in enumerate(CONTEXTS_78)}

BASES = "ACGT"
INT2BASE = np.array(list("ACGTN"))  # supports ref_base==4
fasta = pysam.FastaFile(str(FASTA))

index = defaultdict(list)
row_counter = 0

for fp in sorted(PARQ.glob("chr*.parquet")):
    tbl  = pq.read_table(fp, columns=["chr", "pos", "ref_base"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()                 # 1-based
    refs = INT2BASE[tbl["ref_base"].to_numpy()]  # A/C/G/T/N

    for chrom, pos, ref_base in zip(chrs, poss, refs):
        # Boundary + N guard
        chrom_len = fasta.get_reference_length(chrom)
        pos = int(pos)

        # Need two bases: pos and pos+1 in 1-based coordinates
        if ref_base == "N" or pos < 1 or pos >= chrom_len:
            index[-1].append(row_counter)
            row_counter += 1
            continue

        # Fetch dinucleotide starting at pos (1-based):
        # pysam fetch uses 0-based, end-exclusive => [pos-1, pos+1)
        ref_dinuc = fasta.fetch(chrom, pos - 1, pos + 1).upper()
        if len(ref_dinuc) != 2 or "N" in ref_dinuc:
            index[-1].append(row_counter)
            row_counter += 1
            continue

        # Optional: sanity check that left base matches parquet ref_base
        if ref_dinuc[0] != ref_base:
            index[-1].append(row_counter)
            row_counter += 1
            continue

        matched = False
        for a1 in BASES:
            for a2 in BASES:
                alt = a1 + a2
                if alt == ref_dinuc:
                    continue
                ctx = f"{ref_dinuc}>{alt}"
                ctx_id = CTX2ID.get(ctx)
                if ctx_id is not None:
                    index[ctx_id].append(row_counter)
                    matched = True

        if not matched:
            index[-1].append(row_counter)

        row_counter += 1

    print(f"Processed {fp.name}, total rows so far: {row_counter:,}")

# guarantee all 0..77 keys exist
for i in range(78):
    index[i] = index.get(i, [])

full_index = {k: np.asarray(v, dtype=np.int32) for k, v in index.items()}
with open(OUT, "wb") as f:
    pickle.dump(full_index, f, protocol=4)

print(f"\n✅ Saved {OUT}")
nonempty = sum(1 for k, v in full_index.items() if k != -1 and len(v) > 0)
print(f"Contexts present: {nonempty}/78 (non-empty pools)")
print(f"Fallback (-1) rows: {len(full_index.get(-1, [])):,}")
