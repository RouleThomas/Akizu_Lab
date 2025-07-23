#!/usr/bin/env python3
import pandas as pd
import pysam
import argparse
from Bio.Seq import reverse_complement


def parquet_to_txt(parquet_path, fasta_path, sample_name, output_path):
    df = pd.read_parquet(parquet_path)
    fasta = pysam.FastaFile(fasta_path)

    # Ensure ref and alt are uppercase strings
    df["ref_base"] = df["ref_base"].astype(str).str.upper()
    df["alt_base"] = df["alt_base"].astype(str).str.upper()
    df["simulated_ref"] = df["ref_base"].copy()

    # Apply reverse complement to simulated_ref if revcomp is True
    df.loc[df["revcomp"], "simulated_ref"] = (
        df.loc[df["revcomp"], "simulated_ref"]
        .apply(lambda x: reverse_complement(x))
    )

    # Fetch true genome reference base
    def fetch_genome_ref(row):
        base = fasta.fetch(row.chr, row.pos - 1, row.pos).upper()
        return base

    df["genome_ref"] = df.apply(fetch_genome_ref, axis=1)

    # Validate matches
    match = df["simulated_ref"] == df["genome_ref"]
    n_matches = match.sum()
    n_total = len(df)
    n_mismatches = (~match).sum()

    print(f"✅ {n_matches} / {n_total} matches between genome and simulated ref_base")
    if n_mismatches > 0:
        print(f"❌ {n_mismatches} mismatches found")
        print("Sample mismatches:")
        print(df.loc[~match, ["chr", "pos", "ref_base", "alt_base", "simulated_ref", "genome_ref", "revcomp"]].head(10))

    # Filter only valid mutations
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
    parser.add_argument("--parquet", required=True, help="Input annotated parquet file")
    parser.add_argument("--fasta", required=True, help="Reference genome FASTA file")
    parser.add_argument("--sample-name", required=True, help="Sample name")
    parser.add_argument("--output", required=True, help="Output txt path")
    args = parser.parse_args()

    parquet_to_txt(args.parquet, args.fasta, args.sample_name, args.output)
