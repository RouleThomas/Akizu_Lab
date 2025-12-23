#!/usr/bin/env python3
import pandas as pd
from pathlib import Path
import numpy as np
import pysam

FASTA = "ref/GRCh38.primary_assembly.genome.fa"

def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT", "TGCA"))[::-1]

def main():
    fasta = pysam.FastaFile(FASTA)

    dfs = []
    for p in sorted(Path("parquet").glob("chr*.parquet")):
        df = pd.read_parquet(
            p,
            columns=["chr","pos","strand","transcript_id","codon_pos"]
        )
        dfs.append(df)
    ex = pd.concat(dfs, ignore_index=True)

    # detect codon_pos convention (should be 0/1/2 from your builder)
    cp_vals = set(pd.to_numeric(ex["codon_pos"], errors="coerce").dropna().astype(int).unique())
    if cp_vals.issubset({0,1,2}):
        boundary_val = 2
    elif cp_vals.issubset({1,2,3}):
        boundary_val = 3
    else:
        raise ValueError(f"Unexpected codon_pos values (sample): {sorted(list(cp_vals))[:20]}")

    counts = {"total": {}, "boundary": {}}

    for chrom, sub in ex.groupby("chr", sort=False):
        # lookup codon_pos by (transcript_id, genomic_pos)
        sub_cp = pd.to_numeric(sub["codon_pos"], errors="coerce").astype("Int64")
        keys = list(zip(sub["transcript_id"].astype(str), sub["pos"].astype(int)))
        key2cp = dict(zip(keys, sub_cp.tolist()))

        for r in sub.itertuples(index=False):
            pos = int(r.pos)
            strand = int(r.strand)
            tx = str(r.transcript_id)

            # genomic dinucleotide at (pos, pos+1)
            try:
                pair = fasta.fetch(chrom, pos-1, pos+1).upper()
                if len(pair) != 2 or "N" in pair:
                    continue
            except Exception:
                continue

            ref2 = pair if strand == 1 else revcomp(pair)
            counts["total"][ref2] = counts["total"].get(ref2, 0) + 1

            # strand-aware: coding-start base is pos on + strand, pos+1 on - strand
            start_pos = pos if strand == 1 else (pos + 1)
            cp = key2cp.get((tx, start_pos), None)
            if cp is not None and int(cp) == boundary_val:
                counts["boundary"][ref2] = counts["boundary"].get(ref2, 0) + 1

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
