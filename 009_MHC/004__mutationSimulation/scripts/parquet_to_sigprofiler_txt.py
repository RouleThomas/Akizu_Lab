#!/usr/bin/env python3
import pandas as pd
import pysam
import argparse
from Bio.Seq import reverse_complement


def parquet_to_txt(parquet_path, fasta_path, sample_name, output_path):
    df = pd.read_parquet(parquet_path)
    fasta = pysam.FastaFile(fasta_path)

    # Sanity checks
    if "ref_base_resolved" not in df.columns or "revcomp" not in df.columns:
        raise ValueError("Input parquet must contain 'ref_base_resolved' and 'revcomp' columns")

    # Clean up columns
    df["ref_base_resolved"] = df["ref_base_resolved"].astype(str).str.upper()
    df["alt_base"] = df["alt_base"].astype(str).str.upper()
    df["chr"] = df["chr"].astype(str)

    # Compute simulated_ref with strand correction
    df["simulated_ref"] = df.apply(
        lambda row: reverse_complement(row["ref_base_resolved"]) if row["revcomp"] else row["ref_base_resolved"],
        axis=1
    )

    # Fetch true genome reference base
    def fetch_genome_ref(row):
        try:
            base = fasta.fetch(row.chr, row.pos - 1, row.pos).upper()
            return reverse_complement(base) if row.revcomp else base
        except Exception as e:
            print(f"⚠️ Failed to fetch ref for {row.chr}:{row.pos} — {e}")
            return "N"

    df["genome_ref"] = df.apply(fetch_genome_ref, axis=1)

    # Match validation
    match = df["simulated_ref"] == df["genome_ref"]
    n_matches = match.sum()
    n_total = len(df)
    n_mismatches = (~match).sum()

    print(f"✅ {n_matches} / {n_total} matches between genome and simulated ref_base")
    if n_mismatches > 0:
        print(f"❌ {n_mismatches} mismatches found")
        print("Sample mismatches:")
        print(df.loc[~match, ["chr", "pos", "ref_base_resolved", "alt_base", "simulated_ref", "genome_ref", "revcomp"]].head(10))

    # Keep only validated
    df_valid = df[match].copy()

    # Format for SigProfiler
    df_txt = pd.DataFrame({
        "Project": "Simu",
        "Sample": sample_name,
        "ID": ".",
        "Genome": "GRCh38",
        "mut_type": "SNP",
        "chrom": df_valid["chr"].str.replace("chr", "", regex=False),
        "pos_start": df_valid["pos"],
        "pos_end": df_valid["pos"],
        "ref": df_valid["genome_ref"],
        "alt": df_valid["alt_base"],
        "Type": "SOMATIC"
    })

    df_txt.to_csv(output_path, sep="\t", index=False)
    print(f"✅ Saved {len(df_txt)} validated mutations to {output_path}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True, help="Input annotated .parquet file")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA (indexed)")
    parser.add_argument("--sample-name", required=True, help="Sample name for SigProfiler")
    parser.add_argument("--output", required=True, help="Output path (.txt)")
    args = parser.parse_args()

    parquet_to_txt(args.parquet, args.fasta, args.sample_name, args.output)
