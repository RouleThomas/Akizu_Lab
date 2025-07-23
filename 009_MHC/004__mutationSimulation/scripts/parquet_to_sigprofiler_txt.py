#!/usr/bin/env python3

import argparse
import pandas as pd
import pyarrow.parquet as pq
from Bio import SeqIO
from Bio.Seq import reverse_complement
from pathlib import Path


INT2BASE = ["A", "C", "G", "T"]


def load_genome(fasta_path):
    """Load reference genome into a dictionary."""
    genome = {}
    for record in SeqIO.parse(fasta_path, "fasta"):
        genome[record.id] = str(record.seq).upper()
    return genome


def get_genome_base(genome, chrom, pos):
    """Return base from genome at 1-based position."""
    seq = genome.get(chrom)
    if seq is None or pos - 1 >= len(seq):
        return "N"
    return seq[pos - 1]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    # Load mutation data
    df = pq.read_table(args.parquet).to_pandas()

    # Load genome
    genome = load_genome(args.fasta)

    # Convert integer ref_base to A/C/G/T
    df["simulated_ref"] = [INT2BASE[int(b)] for b in df["ref_base"]]

    # Reverse complement if revcomp is True
    df["simulated_ref"] = [
        reverse_complement(b) if rev else b
        for b, rev in zip(df["simulated_ref"], df["revcomp"])
    ]

    # Extract reference base from genome
    df["genome_ref"] = [
        get_genome_base(genome, chrom, pos)
        for chrom, pos in zip(df["chr"], df["pos"])
    ]

    # Compare genome_ref and simulated_ref
    df["match"] = df["genome_ref"] == df["simulated_ref"]
    n_matches = df["match"].sum()
    n_total = len(df)

    print(f"✅ {n_matches} / {n_total} matches between genome and simulated ref_base")
    print(f"❌ {n_total - n_matches} mismatches found")

    mismatches = df[~df["match"]]
    if not mismatches.empty:
        print(mismatches[["chr", "pos", "ref_base", "alt_base", "simulated_ref", "genome_ref", "revcomp"]].head(10))

    # Save only validated mutations
    df_valid = df[df["match"]].copy()

    # Format for SigProfiler: CHROM POS REF ALT SAMPLE
    df_out = df_valid[["chr", "pos", "simulated_ref", "alt_base"]].copy()
    df_out["sample"] = args.sample_name
    df_out.to_csv(args.output, sep="\t", header=False, index=False)

    print(f"✅ Saved {len(df_out)} validated mutations to {args.output}")


if __name__ == "__main__":
    main()
