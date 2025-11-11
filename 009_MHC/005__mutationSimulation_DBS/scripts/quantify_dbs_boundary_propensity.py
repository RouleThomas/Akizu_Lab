#!/usr/bin/env python3
# quantify_dbs_boundary_propensity.py
import pandas as pd
from pathlib import Path
import numpy as np
import pysam

# Inputs:
# - exome parquet with at least: chr, pos, strand, codon_index (0,1,2 for positions within codon)
#   (same source you used to build indices)
# - reference FASTA
EXOME_PARQUET = "parquet/chr*.parquet"        # glob over your per-chr parquet
FASTA = "ref/GRCh38.primary_assembly.genome.fa"

def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT","TGCA"))[::-1]

def main():
    fasta = pysam.FastaFile(FASTA)
    # load just what we need
    dfs = []
    for p in sorted(Path(".").glob(EXOME_PARQUET)):
        df = pd.read_parquet(p, columns=["chr","pos","strand","codon_index","ref_base"])
        dfs.append(df)
    ex = pd.concat(dfs, ignore_index=True)

    # we’ll count dinucleotides starting at each position (pos, pos+1) on the GENOMIC strand,
    # but report the 5'->3' REF pair on the TRANSCRIPT/coding strand so all codon_index are comparable.
    counts = { "total": {}, "boundary": {} }   # dicts: ref2mer -> count

    # vectorized chunking for speed
    for chrom, sub in ex.groupby("chr"):
        # fetch a big chunk once: make a boolean mask to avoid going out of contig
        # Here we just do row-wise fetch (simpler & fine for a one-off analysis)
        for _, r in sub.iterrows():
            pos = int(r.pos)
            strand = int(r.strand)   # +1 or -1
            codon_idx = int(r.codon_index)  # 0,1,2

            try:
                pair = fasta.fetch(chrom, pos-1, pos+1).upper()  # 1-based pos; fetch [pos, pos+1] -> 0-based [pos-1, pos+1)
                if len(pair) != 2 or "N" in pair: 
                    continue
            except:
                continue

            # orient to coding strand
            ref2 = pair if strand == 1 else revcomp(pair)

            counts["total"][ref2] = counts["total"].get(ref2, 0) + 1
            if codon_idx == 2:   # starting at codon pos 3 -> spans into next codon (two AA)
                counts["boundary"][ref2] = counts["boundary"].get(ref2, 0) + 1

    # assemble table for the 10 possible ref two-mers we actually see in DBS contexts
    wanted = ["AC","AT","CC","CG","CT","GC","TA","TC","TG","TT"]
    rows = []
    for ref2 in wanted:
        tot = counts["total"].get(ref2, 0)
        bnd = counts["boundary"].get(ref2, 0)
        prop = (bnd / tot) if tot else np.nan
        rows.append({"ref2": ref2, "count_total": tot, "count_boundary": bnd, "prop_boundary": prop})

    out = pd.DataFrame(rows).sort_values("prop_boundary", ascending=False)
    print(out.to_string(index=False))
    out.to_csv("dbs_ref2_boundary_propensity.tsv", sep="\t", index=False)
    print("\nSaved: dbs_ref2_boundary_propensity.tsv")

if __name__ == "__main__":
    main()
