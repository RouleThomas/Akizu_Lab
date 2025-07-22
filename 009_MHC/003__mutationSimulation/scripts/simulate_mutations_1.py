#!/usr/bin/env python3
"""
simulate_mutations_1.py
Simulates SBS mutations from a given signature across an exome using a context index,
ensuring exact reference base match by querying the genome immediately.
"""

import argparse, pickle
import numpy as np
import pandas as pd
from pathlib import Path
import pyarrow as pa, pyarrow.parquet as pq
from Bio import SeqIO


def load_signatures(path):
    df = pd.read_csv(path, sep="\t", index_col=0)
    df.columns = [col.upper() for col in df.columns]
    return df


def get_signature_probs(df, signature):
    if signature not in df.columns:
        raise ValueError(f"Signature '{signature}' not found in provided file.")
    return df[signature.upper()].values


class FastaGenome:
    def __init__(self, fasta_path):
        self.ref = {}
        for rec in SeqIO.parse(fasta_path, "fasta"):
            self.ref[rec.id] = str(rec.seq).upper()

    def fetch(self, chrom, pos):
        return self.ref[chrom][pos - 1]  # 1-based to 0-based


class ExomeSampler:
    def __init__(self, exome_dir, context_index_path):
        self.tables = []
        for f in sorted(Path(exome_dir).glob("*.parquet")):
            self.tables.append(pq.read_table(f))
        self.exome = pa.concat_tables(self.tables)

        with open(context_index_path, "rb") as f:
            self.ctx_index = pickle.load(f)

    def sample(self, probs, n_mutations, rng, genome):
        sampled = []
        attempts = 0
        max_attempts = n_mutations * 20

        context_ids = np.arange(len(probs))

        while len(sampled) < n_mutations and attempts < max_attempts:
            ctx = rng.choice(context_ids, p=probs)
            rows = self.ctx_index.get(ctx, [])

            if len(rows) == 0:
                attempts += 1
                continue

            row = int(rows[rng.integers(len(rows))])
            r = self.exome.slice(row, 1).to_pandas().iloc[0]
            chrom, pos, ref_base, alt_base, rev = r.chr, r.pos, r.ref_base, r.alt_base, r.revcomp

            try:
                genome_ref = genome.fetch(chrom, pos)
            except Exception:
                attempts += 1
                continue

            if genome_ref != ref_base:
                attempts += 1
                continue

            sampled.append(row)
            attempts += 1

        print(f"✅ Simulated: {len(sampled)} / {n_mutations} mutations after {attempts} attempts")
        return pa.concat_tables([self.exome.slice(i, 1) for i in sampled])


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--signatures", required=True)
    parser.add_argument("--signature-name", required=True)
    parser.add_argument("--n", type=int, required=True)
    parser.add_argument("--exome-dir", required=True)
    parser.add_argument("--context-index", required=True)
    parser.add_argument("--fasta", required=True)
    parser.add_argument("--out", required=True)
    parser.add_argument("--seed", type=int, default=42)
    args = parser.parse_args()

    sig_df = load_signatures(args.signatures)
    probs = get_signature_probs(sig_df, args.signature_name)

    genome = FastaGenome(args.fasta)
    sampler = ExomeSampler(args.exome_dir, args.context_index)

    rng = np.random.default_rng(args.seed)
    sampled_tbl = sampler.sample(probs, args.n, rng, genome)

    pq.write_table(sampled_tbl, args.out, compression="zstd")
    print(f"✅ Saved: {args.out} ({sampled_tbl.num_rows} mutations)")


if __name__ == "__main__":
    main()
