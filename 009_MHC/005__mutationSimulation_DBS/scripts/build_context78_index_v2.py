#!/usr/bin/env python3
"""
Build DBS context index (0..77 and -1) in the SAME representation
as your signature file (i.e., keep forward-strand ref doublets; NO
pyrimidine-left canonicalization).
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
fasta = pysam.FastaFile(str(FASTA))

index = defaultdict(list)
row_counter = 0

for fp in sorted(PARQ.glob("chr*.parquet")):
    tbl  = pq.read_table(fp, columns=["chr", "pos", "ref_base"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()
    refs = np.take(list(BASES), tbl["ref_base"].to_numpy())  # 'A','C','G','T'

    for chrom, pos, ref in zip(chrs, poss, refs):
        # forward-strand 2-bp reference doublet (NO canonicalization)
        ref_dinuc = fasta.fetch(chrom, pos - 1, pos + 1).upper()
        if len(ref_dinuc) != 2 or "N" in ref_dinuc:
            index[-1].append(row_counter)
            row_counter += 1
            continue

        # attach this row to ALL mutation classes ref_dinuc > ALT (ALT != ref)
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
nonempty = sum(int(len(v) > 0) for v in full_index.values() if k != -1)
print(f"Contexts present: {nonempty}/78 (non-empty pools)")
