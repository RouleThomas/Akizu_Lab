#!/usr/bin/env python3
import pandas as pd
import pysam
import argparse

def parquet_to_txt(parquet, fasta_path, sample_name, output_path):
    df = pd.read_parquet(parquet)
    fasta = pysam.FastaFile(fasta_path)

    # Ensure ref_base is uppercase character, not int (already handled correctly in simulation now)
    df["simulated_ref"] = df["ref_base"].str.upper()

    # Fetch the actual base from the genome FASTA
    df["genome_ref"] = df.apply(lambda r: fasta.fetch(r.chr, r.pos - 1, r.pos).upper(), axis=1)

    # Keep only valid matches
    match = df["genome_ref"] == df["simulated_ref"]
    mismatches = (~match).sum()
    print(f"✅ {match.sum()} / {len(df)} matches")
    if mismatches > 0:
        print(f"❌ {mismatches} mismatches")
        print("Sample mismatches:")
        print(df.loc[~match, ["chr", "pos", "simulated_ref", "genome_ref", "alt_base", "revcomp"]].head())

    df_valid = df[match]

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
