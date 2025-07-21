#!/usr/bin/env python3
import pandas as pd
import pysam
import argparse
from Bio.Seq import reverse_complement


def parquet_to_txt(parquet, fasta_path, sample_name, output_path):
    df = pd.read_parquet(parquet)
    fasta = pysam.FastaFile(fasta_path)

    # Ensure simulated_ref and alt_base are upper-case strings
    df["simulated_ref"] = df["ref_base"].astype(str).str.upper()
    df["alt_base"] = df["alt_base"].astype(str).str.upper()

    def get_genome_ref(row):
        base = fasta.fetch(row.chr, row.pos - 1, row.pos).upper()
        return reverse_complement(base) if row.revcomp else base

    df["genome_ref"] = df.apply(get_genome_ref, axis=1)

    # Match validation
    match = df["genome_ref"] == df["simulated_ref"]
    mismatches = (~match).sum()
    print(f"✅ {match.sum()} / {len(df)} matches")
    if mismatches > 0:
        print(f"❌ {mismatches} mismatches")
        print("Sample mismatches:")
        print(df.loc[~match, ["chr", "pos", "ref_base", "alt_base", "simulated_ref", "genome_ref", "revcomp"]].head(10))

    df_valid = df[match].copy()

    df_txt = pd.DataFrame({
        "Project": "Simu",
        "Sample": sample_name,
        "ID": ".",
        "Genome": "GRCh38",
        "mut_type": "SNP",
        "chrom": df_valid["chr"].str.replace("chr", "", regex=True),
        "pos_start": df_valid["pos"],
        "pos_end": df_valid["pos"],
        "ref": df_valid["genome_ref"],
        "alt": df_valid["alt_base"],
        "Type": "SOMATIC"
    })

    df_txt.to_csv(output_path, sep="\t", index=False)
    print(f"✅ Saved {len(df_txt)} mutations to {output_path}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--output", required=True)
    args = parser.parse_args()

    parquet_to_txt(args.parquet, args.fasta, args.sample_name, args.output)
