#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
import pysam
from Bio.Seq import reverse_complement
from pathlib import Path

# Define COSMIC-style 96 trinucleotide contexts
BASES = "ACGT"
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref in "CT"
    for alt in ("AGT" if ref == "C" else "ACG")
    for l in BASES
    for r in BASES
]
CONTEXT2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}

def get_context(trinuc, ref, alt, revcomp_flag):
    if revcomp_flag:
        trinuc = reverse_complement(trinuc)
        ref = trinuc[1]
        alt = reverse_complement(alt)
    return f"{trinuc[0]}[{ref}>{alt}]{trinuc[2]}"

def main(parquet_path, fasta_path, output_pdf, sample_name):
    df = pd.read_parquet(parquet_path)
    fasta = pysam.FastaFile(fasta_path)

    INT2BASE = np.array(["A", "C", "G", "T"])

    if np.issubdtype(df["ref_base"].dtype, np.integer):
        df["ref_base"] = INT2BASE[df["ref_base"].values]
    if np.issubdtype(df["alt_base"].dtype, np.integer):
        df["alt_base"] = INT2BASE[df["alt_base"].values]

    contexts = []
    mismatches = 0
    mismatch_debug = []

    for _, row in df.iterrows():
        chrom, pos = row["chr"], row["pos"]
        ref, alt = row["ref_base"], row["alt_base"]
        rev = row["revcomp"]

        try:
            trinuc = fasta.fetch(chrom, pos - 2, pos + 1).upper()
            genome_ref = fasta.fetch(chrom, pos - 1, pos).upper()
        except Exception as e:
            contexts.append(None)
            mismatch_debug.append((chrom, pos, "fetch_fail", str(e)))
            mismatches += 1
            continue


        # sanity check
        if ref != genome_ref:
            contexts.append(None)
            mismatch_debug.append((chrom, pos, "ref_mismatch", f"expected: {ref}, got: {genome_ref}"))
            mismatches += 1
            continue

        ctx = get_context(trinuc, ref, alt, revcomp_flag=rev)
        contexts.append(ctx if ctx in CONTEXT2ID else None)

    df["context96"] = contexts
    df_valid = df.dropna(subset=["context96"])

    print(f"✅ {len(df_valid)} valid mutations")
    print(f"❌ {mismatches} mismatches with genome reference")
    if mismatch_debug:
        print("🔎 Sample mismatches (first 10):")
        for i, m in enumerate(mismatch_debug[:10]):
            print(f"  {i+1}. {m}")

    if len(df_valid) == 0:
        print("⚠️ No valid mutations — skipping plot.")
        return

    counts = df_valid["context96"].value_counts().reindex(CONTEXTS_96, fill_value=0)
    percent = 100 * counts / counts.sum()
    bar_colors = [ 
        {"C>A": "#00BFC4", "C>G": "black", "C>T": "red", "T>A": "gray", "T>C": "green", "T>G": "salmon"}[ctx[2:5]]
        for ctx in CONTEXTS_96
    ]

    plt.figure(figsize=(20, 6))
    plt.bar(range(96), percent, color=bar_colors)
    plt.xticks(range(96), CONTEXTS_96, rotation=90, fontsize=5)
    plt.ylabel("Percentage of Single Base Substitutions")
    plt.title(sample_name)
    plt.tight_layout()
    plt.savefig(output_pdf)
    print(f"📊 Plot saved to {output_pdf}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--parquet", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--sample-name", required=True)
    parser.add_argument("--output-pdf", required=True)
    args = parser.parse_args()

    main(args.parquet, args.fasta, args.output_pdf, args.sample_name)
