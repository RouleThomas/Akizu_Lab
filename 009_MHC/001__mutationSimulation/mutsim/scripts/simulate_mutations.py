#!/usr/bin/env python3
# simulate_mutations.py
"""
Fast SBS-signature mutation simulator on a pre-built exome Parquet table.
See header comments for required input artefacts.
"""
from __future__ import annotations

import argparse, json
import pickle
from pathlib import Path
from typing import Dict, Sequence

import numpy as np
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
from Bio.Data import CodonTable

# ---------------------------------------------------------------------
CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
STOP_CODONS = {"TAA", "TAG", "TGA"}
BASES = np.array(["A", "C", "G", "T"], dtype="U1")

CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alts in [("C", "AGT"), ("T", "ACG")]
    for alt in alts
    for l in "ACGT"
    for r in "ACGT"
]
CTX2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}

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
        self.arrow_table = pa.concat_tables([
            pq.read_table(fp) for fp in sorted(Path(parquet_dir).glob("chr*.parquet"))
        ])
        # promote string columns
        new_cols, new_schema = [], []
        for field, col in zip(self.arrow_table.schema, self.arrow_table.columns):
            if pa.types.is_string(field.type):
                col = col.cast(pa.large_string())
                field = pa.field(field.name, pa.large_string())
            new_cols.append(col); new_schema.append(field)
        self.arrow_table = pa.Table.from_arrays(new_cols, schema=pa.schema(new_schema))
        with open(context_index_pkl, "rb") as fh:
            self.context_index: Dict[int, np.ndarray] = pickle.load(fh)

    def sample(self, signature_probs: np.ndarray, n: int, rng=None) -> pa.Table:
        rng    = rng or np.random.default_rng()
        draws  = rng.multinomial(n, signature_probs)
        picked = []
        for ctx_id, k in enumerate(draws):
            if k:
                pool = self.context_index[ctx_id]
                picked.extend(rng.choice(pool, k, replace=False))
        return self.arrow_table.take(picked)

    def annotate(self, tbl: pa.Table, rng=None) -> pd.DataFrame:
        rng = rng or np.random.default_rng()
        pdf = tbl.to_pandas()
        offsets = rng.integers(1, 4, size=len(pdf))
        pdf["alt_base"] = BASES[(pdf["ref_base"] + offsets) % 4]

        mutated, aa_alt, cons = [], [], []
        for ref_codon, idx, alt in zip(pdf.ref_codon, pdf.codon_index, pdf.alt_base):
            cod = list(ref_codon); cod[idx] = alt; new = "".join(cod)
            mutated.append(new)
            if new in STOP_CODONS:
                aa_alt.append("*"); cons.append("stop")
            else:
                aa_new = CODON_TABLE.get(new)
                aa_ref = CODON_TABLE.get(ref_codon)
                if aa_new is None:
                    aa_alt.append(None); cons.append("unknown")
                elif aa_new == aa_ref:
                    aa_alt.append(aa_new); cons.append("synonymous")
                else:
                    aa_alt.append(aa_new); cons.append("missense")

        pdf["mut_codon"]   = mutated
        pdf["alt_aa"]      = aa_alt
        pdf["consequence"] = cons
        return pdf

    @staticmethod
    def summarise(df: pd.DataFrame) -> dict[str,float]:
        return {
            "n_total":      len(df),
            "frac_syn":     float((df.consequence=="synonymous").mean()),
            "frac_stop":    float((df.consequence=="stop").mean()),
            "frac_missense":float((df.consequence=="missense").mean()),
        }

def main():
    p = argparse.ArgumentParser(
        description="Simulate mutations from an SBS signature"
    )
    p.add_argument("--signatures",     required=True)
    p.add_argument("--signature-name", required=True)
    p.add_argument("--n",              type=int, required=True)
    p.add_argument("--exome-dir",      required=True)
    p.add_argument("--context-index",  required=True)
    p.add_argument("--out",            required=True)
    p.add_argument("--seed",           type=int, default=0)
    args = p.parse_args()

    sig_df = load_signatures(args.signatures)
    probs  = get_signature_probs(sig_df, args.signature_name)
    rng    = np.random.default_rng(args.seed)
    samp   = ExomeSampler(args.exome_dir, args.context_index)
    tbl    = samp.sample(probs, args.n, rng)
    df     = samp.annotate(tbl, rng)
    pq.write_table(pa.Table.from_pandas(df), Path(args.out), compression="zstd")
    print(json.dumps(samp.summarise(df)), flush=True)

if __name__ == "__main__":
    main()

