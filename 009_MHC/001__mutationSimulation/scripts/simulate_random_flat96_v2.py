#!/usr/bin/env python3
# simulate_random_flat96.py

import argparse, json, pickle, re
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import pyarrow as pa, pyarrow.parquet as pq
from Bio.Data import CodonTable

from annotate_damage2 import DBNSFP

# Genetic code tables
CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
STOP_CODONS = {"TAA", "TAG", "TGA"}
BASES = np.array(["A", "C", "G", "T"], dtype="U1")

# 96 trinucleotide mutation contexts
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alts in [("C", "AGT"), ("T", "ACG")]
    for alt in alts
    for l in "ACGT"
    for r in "ACGT"
]
CTX2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}

def get_flat_uniform_probs() -> np.ndarray:
    return np.full(96, 1 / 96)

class ExomeSampler:
    def __init__(self, parquet_dir: str | Path, context_index_pkl: str | Path):
        self.arrow_table = pa.concat_tables([
            pq.read_table(fp) for fp in sorted(Path(parquet_dir).glob("chr*.parquet"))
        ])
        # Convert strings to large_string
        new_cols, new_schema = [], []
        for field, col in zip(self.arrow_table.schema, self.arrow_table.columns):
            if pa.types.is_string(field.type):
                col = col.cast(pa.large_string())
                field = pa.field(field.name, pa.large_string())
            new_cols.append(col)
            new_schema.append(field)
        self.arrow_table = pa.Table.from_arrays(new_cols, schema=pa.schema(new_schema))

        with open(context_index_pkl, "rb") as fh:
            self.context_index: Dict[int, np.ndarray] = pickle.load(fh)

    def sample(self, probs: np.ndarray, n: int, rng=None) -> tuple[pa.Table, np.ndarray]:
        rng = rng or np.random.default_rng()
        draws = rng.multinomial(n, probs)
        picked, picked_ctx_ids = [], []
        for ctx_id, k in enumerate(draws):
            if k:
                pool = self.context_index[ctx_id]
                sampled_rows = rng.choice(pool, k, replace=False)
                picked.extend(sampled_rows)
                picked_ctx_ids.extend([ctx_id] * k)
        return self.arrow_table.take(picked), np.array(picked_ctx_ids)

    def annotate(self, tbl: pa.Table, ctx_ids: np.ndarray) -> pd.DataFrame:
        pdf = tbl.to_pandas()
        pdf["context_id"] = ctx_ids
        alt_bases = [re.search(r"\[(.?)>(.?)\]", CONTEXTS_96[ctx_id]).group(2) for ctx_id in ctx_ids]
        pdf["alt_base"] = alt_bases

        mutated, aa_alt, cons = [], [], []
        for ref_codon, idx, alt in zip(pdf.ref_codon, pdf.codon_index, pdf["alt_base"]):
            cod = list(ref_codon)
            cod[idx] = alt
            new = "".join(cod)
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

def summarize_output(df: pd.DataFrame) -> dict:
    consequence_counts = df["consequence"].value_counts().to_dict()
    scores = {
        "sift4g_score": df["sift4g_score"].dropna().tolist(),
        "polyphen2_hdiv_score": df["polyphen2_hdiv_score"].dropna().tolist(),
        "cadd_phred": df["cadd_phred"].dropna().tolist(),
    }
    return {
        "n_mutations": len(df),
        "consequences": consequence_counts,
        "score_counts": {k: len(v) for k, v in scores.items()},
        "score_means": {k: round(np.mean(v), 3) if v else None for k, v in scores.items()},
    }

def main():
    p = argparse.ArgumentParser(description="Simulate mutations with uniform flat 1/96 distribution across trinucleotide contexts")
    p.add_argument("--n",              type=int, required=True)
    p.add_argument("--exome-dir",      required=True)
    p.add_argument("--context-index",  required=True)
    p.add_argument("--rep",            type=int, required=True)
    p.add_argument("--out",            required=True)
    p.add_argument("--seed",           type=int, default=0)
    p.add_argument("--dbnsfp",         required=True)
    args = p.parse_args()

    rng = np.random.default_rng(args.seed)
    probs = get_flat_uniform_probs()

    # Sample and annotate
    sampler = ExomeSampler(args.exome_dir, args.context_index)
    tbl, ctx_ids = sampler.sample(probs, args.n, rng)
    df = sampler.annotate(tbl, ctx_ids)

    # Add SIFT/PolyPhen/CADD
    db = DBNSFP(Path(args.dbnsfp))
    annotations = {k: [] for k in [
        "polyphen2_hdiv_score", "polyphen2_hdiv_pred",
        "sift4g_score", "sift4g_pred", "cadd_phred"]}

    for chrom, pos, r_idx, a in zip(df.chr, df.pos, df.ref_base, df.alt_base):
        ref = "ACGT"[int(r_idx)]
        result = db.query(str(chrom), int(pos), ref, str(a))
        for key in annotations:
            annotations[key].append(result[key])
    for k, v in annotations.items():
        df[k] = v

    # Save output
    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    pq.write_table(pa.Table.from_pandas(df), out_path, compression="zstd")

    json_path = out_path.with_suffix(".summary.json")
    with open(json_path, "w") as f:
        json.dump(summarize_output(df), f, indent=2)

    print(f"✅ Saved: {out_path} and {json_path}")

if __name__ == "__main__":
    main()
