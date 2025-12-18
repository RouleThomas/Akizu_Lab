#!/usr/bin/env python3
from __future__ import annotations

import argparse, json, pickle, re
from pathlib import Path
from typing import Dict

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq
from Bio.Data import CodonTable

CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
STOP_CODONS = {"TAA", "TAG", "TGA"}

# IMPORTANT: include N for safety (ref_base can be 4)
INT2BASE = ["A", "C", "G", "T", "N"]

# COSMIC SBS96 contexts in canonical C/T representation
CONTEXTS_96 = [
    f"{l}[{ref}>{alt}]{r}"
    for ref, alts in [("C", "AGT"), ("T", "ACG")]
    for alt in alts
    for l in "ACGT"
    for r in "ACGT"
]
CTX2ID = {ctx: i for i, ctx in enumerate(CONTEXTS_96)}

# Fast revcomp for single bases / short strings
_RC = str.maketrans("ACGTNacgtn", "TGCANtgcan")
def revcomp(seq: str) -> str:
    return seq.translate(_RC)[::-1]

# Regex to parse contexts like A[C>T]G
CTX_RE = re.compile(r"^([ACGT])\[([CT])>([ACGT])\]([ACGT])$")


def load_signatures(path: str | Path) -> pd.DataFrame:
    df = pd.read_csv(path, sep="\t", index_col=0)
    # enforce SBS96 order
    df = df.reindex(CONTEXTS_96)
    # normalize each signature column to sum to 1
    return df / df.sum(axis=0)

def get_signature_probs(sig_df: pd.DataFrame, signature: str) -> np.ndarray:
    if signature not in sig_df.columns:
        raise KeyError(f"Signature {signature} not in file")
    vec = sig_df[signature].to_numpy(dtype="float64")
    s = float(vec.sum())
    if s <= 0:
        raise ValueError(f"Signature {signature} has zero/invalid probabilities")
    return vec / s


class ExomeSampler:
    def __init__(self, parquet_dir: str | Path, context_index_pkl: str | Path):
        parquet_dir = Path(parquet_dir)

        # Keep row order consistent with context index builder: sorted chr*.parquet
        tables = [pq.read_table(fp) for fp in sorted(parquet_dir.glob("chr*.parquet"))]
        if not tables:
            raise FileNotFoundError(f"No chr*.parquet found in {parquet_dir}")
        self.arrow_table = pa.concat_tables(tables)

        # Ensure string columns are large_string to avoid pyarrow/pandas edge issues
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

        # drop fallback bucket if present
        if -1 in self.context_index:
            del self.context_index[-1]

    def sample(self, signature_probs: np.ndarray, n: int, rng=None) -> tuple[pa.Table, np.ndarray]:
        rng = rng or np.random.default_rng()

        draws = rng.multinomial(n, signature_probs)
        picked = []
        picked_ctx_ids = []

        for ctx_id, k in enumerate(draws):
            if k <= 0:
                continue
            pool = self.context_index.get(ctx_id)
            if pool is None or len(pool) == 0:
                continue

            # If k > pool size, you must sample with replacement, otherwise it crashes.
            replace = k > len(pool)
            sampled_rows = rng.choice(pool, k, replace=replace)
            picked.extend(sampled_rows.tolist() if hasattr(sampled_rows, "tolist") else list(sampled_rows))
            picked_ctx_ids.extend([ctx_id] * k)

        return self.arrow_table.take(pa.array(picked, type=pa.int64())), np.array(picked_ctx_ids, dtype=np.int32)

    def annotate(self, tbl: pa.Table, ctx_ids: np.ndarray) -> pd.DataFrame:
        pdf = tbl.to_pandas()
        pdf["context_id"] = ctx_ids

        # Parse contexts -> canonical ref (C/T) and canonical alt (A/C/G/T)
        ref_from_ctx = np.empty(len(ctx_ids), dtype="U1")
        alt_from_ctx = np.empty(len(ctx_ids), dtype="U1")

        for i, ctx_id in enumerate(ctx_ids):
            ctx = CONTEXTS_96[int(ctx_id)]
            m = CTX_RE.match(ctx)
            if not m:
                ref_from_ctx[i] = "N"
                alt_from_ctx[i] = "N"
                continue
            _, ref, alt, _ = m.groups()
            ref_from_ctx[i] = ref
            alt_from_ctx[i] = alt

        pdf["ref_from_ctx"] = ref_from_ctx
        pdf["alt_base"] = alt_from_ctx

        # Resolve ref_base integer -> base
        rb = pdf["ref_base"].astype(int).to_numpy()
        rb = np.clip(rb, 0, 4)  # protect against weird values
        ref_base_resolved = np.array([INT2BASE[x] for x in rb], dtype="U1")
        pdf["ref_base_resolved"] = ref_base_resolved

        # Determine whether the genomic ref is opposite to canonical (C/T) context center:
        # If canonical ref is C/T but actual genomic base is A/G, then we are on the reverse complement.
        # In practice: if ref_from_ctx != actual_ref_base then we need to revcomp the alt base.
        rev = (pdf["ref_from_ctx"].to_numpy(dtype="U1") != ref_base_resolved)
        pdf["revcomp"] = rev

        # Flip alt base if revcomp, and compute simulated_ref consistently
        alt = pdf["alt_base"].to_numpy(dtype="U1")
        alt_flipped = np.array([revcomp(a) if r else a for a, r in zip(alt, rev)], dtype="U1")
        pdf["alt_base"] = alt_flipped

        sim_ref = np.array([revcomp(b) if r else b for b, r in zip(ref_base_resolved, rev)], dtype="U1")
        pdf["simulated_ref"] = sim_ref

        # --- Codon mutation / consequence ---
        # Use codon_pos (0/1/2) if present; else fall back to old codon_index-as-pos behavior.
        if "codon_pos" in pdf.columns:
            codon_pos = pdf["codon_pos"].astype(int).to_numpy()
        else:
            codon_pos = pdf["codon_index"].astype(int).to_numpy()  # backward compatibility

        # guard: force 0/1/2
        codon_pos = np.mod(codon_pos, 3)

        mutated, aa_alt, cons = [], [], []
        for ref_codon, pos_in_codon, alt_base in zip(pdf["ref_codon"], codon_pos, pdf["alt_base"]):
            ref_codon = str(ref_codon).upper()
            if len(ref_codon) != 3 or "N" in ref_codon:
                mutated.append(ref_codon)
                aa_alt.append(None)
                cons.append("unknown")
                continue

            cod = list(ref_codon)
            cod[int(pos_in_codon)] = str(alt_base)
            new = "".join(cod)
            mutated.append(new)

            aa_ref = CODON_TABLE.get(ref_codon, "*" if ref_codon in STOP_CODONS else None)
            aa_new = CODON_TABLE.get(new, "*" if new in STOP_CODONS else None)

            if aa_new is None or aa_ref is None:
                aa_alt.append(None)
                cons.append("unknown")
            elif aa_new == "*":
                aa_alt.append("*")
                cons.append("stop")
            elif aa_new == aa_ref:
                aa_alt.append(aa_new)
                cons.append("synonymous")
            else:
                aa_alt.append(aa_new)
                cons.append("missense")

        pdf["mut_codon"] = mutated
        pdf["alt_aa"] = aa_alt
        pdf["consequence"] = cons

        return pdf

    @staticmethod
    def summarise(df: pd.DataFrame) -> dict[str, float]:
        return {
            "n_total": len(df),
            "frac_syn": float((df.consequence == "synonymous").mean()) if len(df) else 0.0,
            "frac_stop": float((df.consequence == "stop").mean()) if len(df) else 0.0,
            "frac_missense": float((df.consequence == "missense").mean()) if len(df) else 0.0,
        }


def main():
    p = argparse.ArgumentParser(description="Simulate mutations from an SBS signature", allow_abbrev=False)
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

    samp = ExomeSampler(args.exome_dir, args.context_index)
    tbl, ctx_ids = samp.sample(probs, args.n, rng)
    df = samp.annotate(tbl, ctx_ids)

    pq.write_table(pa.Table.from_pandas(df), Path(args.out), compression="zstd")
    print(json.dumps(samp.summarise(df)), flush=True)

if __name__ == "__main__":
    main()
