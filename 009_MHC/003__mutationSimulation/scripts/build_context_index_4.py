#!/usr/bin/env python3
import pickle
from collections import defaultdict
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
import pysam
from Bio.Seq import reverse_complement

HOME = Path.cwd()
PARQ_DIR = HOME / "parquet"
FASTA = HOME / "ref" / "GRCh38.primary_assembly.genome.fa"
OUT = HOME / "indices" / "context96.pkl"
OUT.parent.mkdir(parents=True, exist_ok=True)

SUBS = [("C", "A"), ("C", "G"), ("C", "T"), ("T", "A"), ("T", "C"), ("T", "G")]
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alt in SUBS
    for l in "ACGT"
    for r in "ACGT"
]
CTX2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}
INT2BASE = np.array(list("ACGT"))

fasta = pysam.FastaFile(str(FASTA))

context_rows = defaultdict(list)
row_counter = 0

for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    tbl = pq.read_table(parquet, columns=["chr", "pos", "ref_base", "strand"])
    chrs, poss = tbl["chr"].to_pylist(), tbl["pos"].to_numpy()
    refs, strands = INT2BASE[tbl["ref_base"].to_numpy()], tbl["strand"].to_numpy()

    for chrom, pos, ref, strand in zip(chrs, poss, refs, strands):
        try:
            tri = fasta.fetch(chrom, pos - 1, pos + 2).upper()
            if len(tri) != 3 or "N" in tri:
                row_counter += 1
                continue

            original_ref = tri[1]
            normalized = original_ref in "AG"
            if normalized:
                tri = reverse_complement(tri)
                ref_base = tri[1]
            else:
                ref_base = original_ref

            if ref_base not in "CT":
                row_counter += 1
                continue

            for alt in "AGT" if ref_base == "C" else "ACG":
                if alt == ref_base:
                    continue
                ctx = f"{tri[0]}[{ref_base}>{alt}]{tri[2]}"
                context_rows[CTX2ID[ctx]].append((row_counter, normalized))

            row_counter += 1

        except Exception:
            row_counter += 1
            continue

    print(f"Processed {parquet.name}, cumulative rows = {row_counter:,}")

full_index = {
    ctx_id: np.array(rows, dtype=[('row', 'i4'), ('revcomp', 'b')])
    for ctx_id, rows in context_rows.items()
}

with open(OUT, "wb") as fh:
    pickle.dump(full_index, fh, protocol=4)

print(f"✅ Saved {OUT} with {len(full_index)} contexts.")
