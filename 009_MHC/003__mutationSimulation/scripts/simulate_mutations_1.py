#!/usr/bin/env python3
import argparse, json, pickle
from pathlib import Path
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
from Bio.Seq import reverse_complement
from Bio.Data import CodonTable
import pysam
import pyarrow.compute as pc



CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
STOP_CODONS = {"TAA", "TAG", "TGA"}
BASES = ["A", "C", "G", "T"]
CONTEXTS_96 = [f"{l}[{ref}>{alt}]{r}" for ref, alts in [("C", "AGT"), ("T", "ACG")] for alt in alts for l in BASES for r in BASES]

def get_context(tri: str) -> tuple[str, str, str, bool]:
    """Normalize context to pyrimidine-centered COSMIC style."""
    if tri[1] in "CT":
        return tri[0], tri[1], tri[2], False
    else:
        rc = reverse_complement(tri)
        return rc[0], rc[1], rc[2], True

def validate_ref(chrom, pos, ref_base, fasta, revcomp) -> bool:
    genome_base = fasta.fetch(chrom, pos - 1, pos).upper()
    if revcomp:
        genome_base = reverse_complement(genome_base)
    return genome_base == ref_base

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--signatures", required=True)
    parser.add_argument("--signature-name", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--exome-dir", required=True)
    parser.add_argument("--context-index", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()

    # Load signature
    sig_df = pd.read_csv(args.signatures, sep="\t", index_col=0)
    sig_vec = sig_df[args.signature_name].reindex(CONTEXTS_96).fillna(0).to_numpy()
    probs = sig_vec / sig_vec.sum()
    rng = np.random.default_rng(args.seed)

    # Load index
    with open(args.context_index, "rb") as fh:
        ctx_index = pickle.load(fh)

    # Load exome
    tables = [pq.read_table(f) for f in sorted(Path(args.exome_dir).glob("chr*.parquet"))]
    exome = pa.concat_tables(tables)
    tbl_df = exome.to_pandas()
    fasta = pysam.FastaFile(args.fasta)

    # Sample mutations
    rows, ctx_ids, revs = [], [], []
    while len(rows) < args.n:
        ctx_id = rng.choice(len(CONTEXTS_96), p=probs)
        entries = ctx_index.get(ctx_id, [])
        if len(entries) == 0:
            continue
        entry = rng.choice(entries)
        row_id, revcomp = entry["row"], entry["revcomp"]
        row = tbl_df.iloc[row_id]
        chrom, pos = row["chr"], row["pos"]
        ctx = CONTEXTS_96[ctx_id]
        ref_base = reverse_complement(ctx[2]) if revcomp else ctx[2]

        if validate_ref(chrom, pos, ref_base, fasta, revcomp):
            rows.append(row_id)
            ctx_ids.append(ctx_id)
            revs.append(revcomp)

    sampled_tbl = exome.take(pa.array(list(map(int, rows)), type=pa.int32()))
    df = sampled_tbl.to_pandas()
    df["context_id"] = ctx_ids
    df["revcomp"] = revs
    df["ref_base"] = [reverse_complement(CONTEXTS_96[i][2]) if r else CONTEXTS_96[i][2] for i, r in zip(ctx_ids, revs)]
    df["alt_base"] = [reverse_complement(CONTEXTS_96[i][4]) if r else CONTEXTS_96[i][4] for i, r in zip(ctx_ids, revs)]

    # Codon and consequence
    codon_col, alt_codon_col, aa_col, consequence_col = [], [], [], []
    for i, row in df.iterrows():
        codon = list(row["ref_codon"])
        idx = row["codon_index"]
        alt = df.at[i, "alt_base"]
        rev = row["revcomp"]
        if rev:
            codon_rc = list(reverse_complement("".join(codon)))
            codon_rc[2 - idx] = alt
            mut = reverse_complement("".join(codon_rc))
        else:
            codon[idx] = alt
            mut = "".join(codon)

        alt_codon_col.append(mut)
        aa_ref = CODON_TABLE.get(row["ref_codon"])
        aa_mut = CODON_TABLE.get(mut)
        if mut in STOP_CODONS:
            consequence_col.append("stop")
            aa_col.append("*")
        elif aa_ref == aa_mut:
            consequence_col.append("synonymous")
            aa_col.append(aa_mut or "")
        elif aa_mut:
            consequence_col.append("missense")
            aa_col.append(aa_mut)
        else:
            consequence_col.append("unknown")
            aa_col.append("")

    df["mut_codon"] = alt_codon_col
    df["alt_aa"] = aa_col
    df["consequence"] = consequence_col


    table = pa.Table.from_pandas(df)
    schema = pa.schema([
        (field.name, pa.large_string() if pa.types.is_string(field.type) else field.type)
        for field in table.schema
    ])
    table = pa.Table.from_pandas(df, schema=schema)
    pq.write_table(table, args.out, compression="zstd")

    print(json.dumps({"n_mutations": len(df)}, indent=2))

if __name__ == "__main__":
    main()
