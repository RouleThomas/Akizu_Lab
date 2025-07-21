#!/usr/bin/env python3
import pickle
from collections import defaultdict
from pathlib import Path
import numpy as np
import pyarrow.parquet as pq
import pysam
from Bio.Seq import reverse_complement

# Settings
HOME = Path.cwd()
PARQ_DIR = HOME / "parquet"
FASTA = HOME / "ref" / "GRCh38.primary_assembly.genome.fa"
OUT = HOME / "indices" / "context96.pkl"
OUT.parent.mkdir(parents=True, exist_ok=True)

# COSMIC-style trinucleotide contexts
SUBS = [("C", "A"), ("C", "G"), ("C", "T"), ("T", "A"), ("T", "C"), ("T", "G")]
CONTEXTS_96 = [f"{l}[{ref}>{alt}]{r}" for ref, alt in SUBS for l in "ACGT" for r in "ACGT"]
CTX2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}
INT2BASE = np.array(list("ACGT"))
fasta = pysam.FastaFile(str(FASTA))

context_rows = defaultdict(list)
row_counter = 0
used_rows = 0

for parquet in sorted(PARQ_DIR.glob("chr*.parquet")):
    tbl = pq.read_table(parquet, columns=["chr", "pos", "ref_base"])
    chrs = tbl["chr"].to_pylist()
    poss = tbl["pos"].to_numpy()
    ref_bases = INT2BASE[tbl["ref_base"].to_numpy()]

    for chrom, pos, ref in zip(chrs, poss, ref_bases):
        try:
            tri = fasta.fetch(chrom, pos - 1, pos + 2).upper()
            if len(tri) != 3 or "N" in tri:
                row_counter += 1
                continue

            center = tri[1]
            if center not in "ACGT":
                row_counter += 1
                continue

            # Normalize to pyrimidine-centered view (COSMIC convention)
            if center in "AG":
                norm_tri = reverse_complement(tri)
                ref_base = norm_tri[1]
                l, r = norm_tri[0], norm_tri[2]
                rev = True
            else:
                ref_base = center
                l, r = tri[0], tri[2]
                rev = False

            # Only use valid ref bases
            if ref_base not in "CT":
                row_counter += 1
                continue

            # Add all possible substitutions for that ref base
            alts = "AGT" if ref_base == "C" else "ACG"
            for alt in alts:
                if alt == ref_base:
                    continue
                ctx = f"{l}[{ref_base}>{alt}]{r}"
                if ctx in CTX2ID:
                    context_rows[CTX2ID[ctx]].append((row_counter, rev))
                    used_rows += 1

            row_counter += 1
        except Exception:
            row_counter += 1
            continue

    print(f"✔️ Processed {parquet.name}, cumulative rows = {row_counter:,}, used = {used_rows:,}")

# Save final context index
final = {
    ctx_id: np.array(rows, dtype=[('row', 'i4'), ('revcomp', 'b')])
    for ctx_id, rows in context_rows.items()
}

with open(OUT, "wb") as f:
    pickle.dump(final, f, protocol=4)

print(f"\n✅ Saved context index: {OUT}")
print(f"✅ Indexed contexts: {len(final)}")
print(f"📌 Total rows used for context simulation: {used_rows:,} out of {row_counter:,}")
