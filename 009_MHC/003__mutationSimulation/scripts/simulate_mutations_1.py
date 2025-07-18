#!/usr/bin/env python3

import argparse, json, pickle
from pathlib import Path
from typing import Dict
import numpy as np, pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
from Bio.Data import CodonTable
from Bio.Seq import reverse_complement

CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
STOP_CODONS = {"TAA", "TAG", "TGA"}

BASES = ["A", "C", "G", "T"]
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alts in [("C", "AGT"), ("T", "ACG")]
    for alt in alts
    for l in BASES
    for r in BASES
]

def load_signatures(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    df = df.reindex(CONTEXTS_96)
    return df / df.sum(axis=0)

def get_signature_probs(sig_df: pd.DataFrame, signature: str) -> np.ndarray:
    if signature not in sig_df.columns:
        raise KeyError(f"Signature {signature} not in file")
    vec = sig_df[signature].to_numpy(dtype="float64")
    return vec / vec.sum()

class ExomeSampler:
    def __init__(self, parquet_dir: str | Path, context_index_pkl: str | Path):
        tables = []
        for fp in sorted(Path(parquet_dir).glob("chr*.parquet")):
            tbl = pq.read_table(fp)
            tbl = tbl.cast(
                pa.schema([
                    pa.field(col.name, pa.large_string() if pa.types.is_string(col.type) else col.type)
                    for col in tbl.schema
                ])
            )
            tables.append(tbl)
        self.arrow_table = pa.concat_tables(tables)

        with open(context_index_pkl, "rb") as fh:
            self.context_index: Dict[int, np.ndarray] = pickle.load(fh)

    def sample(self, signature_probs: np.ndarray, n: int, rng=None):
        rng = rng or np.random.default_rng()
        draws = rng.multinomial(n, signature_probs)
        picked_rows, picked_ctx_ids, revcomp_flags = [], [], []

        for ctx_id, k in enumerate(draws):
            if k:
                sampled = rng.choice(self.context_index[ctx_id], k, replace=False)
                picked_rows.extend([s["row"] for s in sampled])
                revcomp_flags.extend([s["revcomp"] for s in sampled])
                picked_ctx_ids.extend([ctx_id] * k)

        return (
            self.arrow_table.take(pa.array(picked_rows, type=pa.int32())),
            np.array(picked_ctx_ids),
            np.array(revcomp_flags, dtype=bool)
        )





    def annotate(self, tbl: pa.Table, ctx_ids: np.ndarray, revcomp_flags: np.ndarray) -> pd.DataFrame:
        df = tbl.to_pandas()
        df["context_id"] = ctx_ids
        df["revcomp"] = revcomp_flags

        mutated_codons = []
        aa_alt = []
        consequences = []
        ref_bases = []
        alt_bases = []

        for i in range(len(df)):
            ref_codon = df.at[i, "ref_codon"]
            codon_index = df.at[i, "codon_index"]
            revcomp = revcomp_flags[i]
            ctx_id = ctx_ids[i]
            context = CONTEXTS_96[ctx_id]

            codon = list(ref_codon)

            if revcomp:
                ref_base = reverse_complement(context[2])
                alt_base = reverse_complement(context[6])
                codon_rc = list(reverse_complement("".join(codon)))
                codon_rc[2 - codon_index] = alt_base
                mutated_codon = reverse_complement("".join(codon_rc))
            else:
                ref_base = context[2]
                alt_base = context[3]
                codon[codon_index] = alt_base
                mutated_codon = "".join(codon)

            ref_bases.append(ref_base)
            alt_bases.append(alt_base)
            mutated_codons.append(mutated_codon)

            if mutated_codon in STOP_CODONS:
                aa_alt.append("*")
                consequences.append("stop")
            else:
                aa_new = CODON_TABLE.get(mutated_codon)
                aa_ref = CODON_TABLE.get(ref_codon)
                if aa_new is None:
                    aa_alt.append(None)
                    consequences.append("unknown")
                elif aa_new == aa_ref:
                    aa_alt.append(aa_new)
                    consequences.append("synonymous")
                else:
                    aa_alt.append(aa_new)
                    consequences.append("missense")

        df["mut_codon"] = mutated_codons
        df["alt_aa"] = aa_alt
        df["consequence"] = consequences
        df["ref_base"] = ref_bases
        df["alt_base"] = alt_bases

        return df





def main():
    p = argparse.ArgumentParser(description="Simulate mutations from an SBS signature")
    p.add_argument("--signatures", required=True)
    p.add_argument("--signature-name", required=True)
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--exome-dir", required=True)
    p.add_argument("--context-index", required=True)
    p.add_argument("--out", required=True)
    p.add_argument("--seed", type=int, default=0)
    args = p.parse_args()

    sig_df = load_signatures(args.signatures)
    probs = get_signature_probs(sig_df, args.signature_name)
    rng = np.random.default_rng(args.seed)

    sampler = ExomeSampler(args.exome_dir, args.context_index)
    tbl, ctx_ids, rev_flags = sampler.sample(probs, args.n, rng)
    df = sampler.annotate(tbl, ctx_ids, rev_flags)

    pq.write_table(pa.Table.from_pandas(df), Path(args.out), compression="zstd")
    print(json.dumps({"n_total": len(df)}, indent=2))

if __name__ == "__main__":
    main()
