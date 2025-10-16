#!/usr/bin/env python3
"""
simulate_mutations_DBS.py — DBS / simulation-only helpers.

- load_signatures(path, context_order_path=None)
- get_signature_probs(sig_df, signature)
- ExomeSampler(exome_dir, context_index_pkl, context_names=None).sample(...)

This module is robust to:
- signature files with BOM/whitespace,
- context index keys being strings (e.g., "AC>CA") or ints (0..N-1 or 1..N),
- sparse/missing contexts in the index (we renormalize over the keys that exist).
"""

from __future__ import annotations
from pathlib import Path
from typing import Dict, Tuple, List, Any
import io, codecs, json, pickle

import numpy as np
import pandas as pd
import pyarrow as pa
import pyarrow.parquet as pq


# ---------------------------- signatures ------------------------------------ #

def load_signatures(path: str | Path, context_order_path: str | Path | None = None) -> pd.DataFrame:
    """Load DBS signature matrix (TSV); first column 'type'/'Type' becomes index.
    If context_order_path is provided, rows are reindexed to that order.
    Columns are normalized to sum to 1."""
    path = Path(path)
    with open(path, "rb") as fh:
        raw = fh.read()
    text = codecs.decode(raw, "utf-8-sig")

    df = pd.read_csv(io.StringIO(text), sep="\t", dtype=str)
    df.columns = [c.strip().replace("\u00a0", " ") for c in df.columns]

    if "type" in df.columns:
        ctx_col = "type"
    elif "Type" in df.columns:
        ctx_col = "Type"
    else:
        raise ValueError(f"{path} must have a first column named 'type' or 'Type'.")

    df[ctx_col] = df[ctx_col].astype(str).str.strip().str.replace("\u00a0", " ", regex=False)
    df = df.set_index(ctx_col)

    df = df.apply(pd.to_numeric, errors="coerce").fillna(0.0)
    col_sums = df.sum(axis=0).replace(0, np.nan)
    df = df.divide(col_sums, axis=1).fillna(0.0)

    if context_order_path is not None:
        with open(context_order_path, "r") as fh:
            order = [ln.strip() for ln in fh if ln.strip()]
        missing = [c for c in order if c not in df.index]
        if missing:
            raise ValueError(
                f"{len(missing)} contexts from {context_order_path} not found in {path}, e.g. {missing[:5]}"
            )
        df = df.reindex(order)

    return df


def get_signature_probs(sig_df: pd.DataFrame, signature: str) -> np.ndarray:
    sig_norm = signature.strip().replace("\u00a0", " ")
    cols_norm = {c.strip().replace("\u00a0", " "): c for c in sig_df.columns}
    if sig_norm not in cols_norm:
        raise KeyError(f"Signature '{signature}' not found. First 12 columns: {list(cols_norm)[:12]}")
    col = cols_norm[sig_norm]
    vec = sig_df[col].to_numpy(dtype="float64")
    s = vec.sum()
    if s <= 0:
        raise ValueError(f"Signature '{signature}' has zero total probability.")
    return vec / s



def load_int_keymap(path):
    m = {}
    with open(path) as f:
        for ln in f:
            ln = ln.strip()
            if not ln: continue
            k, v = ln.split(None, 1)
            m[int(k)] = v.strip()
    return m


# ---------------------------- sampler --------------------------------------- #

def _norm_name(x: Any) -> str:
    return str(x).strip().replace("\u00a0", " ")

class ExomeSampler:
    """Sample DBS contexts from an exome table using a context index.

    context_index_pkl can map:
      { int_id -> np.ndarray[row_idx] } or { "AC>CA" -> np.ndarray[row_idx] }.

    If the index uses string keys, pass `context_names` in the exact order used
    to build the index (typically signatures/context_signature_list_DBS.txt).
    """

    def __init__(
        self,
        parquet_dir: str | Path,
        context_index_pkl: str | Path,
        context_names: List[str] | None = None,
    ):
        parquet_dir = Path(parquet_dir)
        shards = sorted(parquet_dir.glob("*.parquet"))
        if not shards:
            raise FileNotFoundError(f"No parquet shards found in {parquet_dir}")

        tables = [pq.read_table(fp) for fp in shards]
        tbl = pa.concat_tables(tables, promote=True)  # warning about promote is harmless

        # ensure large_string for safety
        new_cols, new_schema = [], []
        for field, col in zip(tbl.schema, tbl.columns):
            if pa.types.is_string(field.type):
                col = col.cast(pa.large_string())
                field = pa.field(field.name, pa.large_string())
            new_cols.append(col)
            new_schema.append(field)
        self.arrow_table = pa.Table.from_arrays(new_cols, schema=pa.schema(new_schema))

        with open(context_index_pkl, "rb") as fh:
            raw_index = pickle.load(fh)

        if not isinstance(raw_index, dict) or not raw_index:
            raise ValueError("Context index is empty or not a dict.")

        if -1 in raw_index:  # drop sentinel if present
            del raw_index[-1]

        # Normalize index to: Dict[key_any, np.ndarray[int]]
        self.raw_index: Dict[Any, np.ndarray] = {k: np.asarray(v, dtype=np.int64) for k, v in raw_index.items()}
        self.context_names = [_norm_name(x) for x in (context_names or [])]
        self._n_rows = self.arrow_table.num_rows

    def _key_to_prob_index(self, key: Any, probs_len: int) -> int | None:
        """Map an index key to the position in the probs vector.
        Supports:
          - str key -> position via self.context_names
          - int key in [0..N-1]
          - int key in [1..N]  (auto shift)
        Returns the position or None if it can't be mapped.
        """
        if isinstance(key, str):
            if not self.context_names:
                return None
            nm = _norm_name(key)
            try:
                return self.context_names.index(nm)
            except ValueError:
                return None
        # treat numpy ints same as python int
        try:
            ki = int(key)
        except Exception:
            return None
        if 0 <= ki < probs_len:
            return ki
        if 1 <= ki <= probs_len:
            return ki - 1  # auto-shift 1-based to 0-based
        return None

    def dump_debug(self, out_json: str | Path, probs_len: int):
        # summary of pools & mapping coverage
        nonempty = sum(int(len(v) > 0) for v in self.raw_index.values())
        # mapping coverage
        ok_map = 0
        for k in self.raw_index.keys():
            if self._key_to_prob_index(k, probs_len) is not None:
                ok_map += 1
        stats = {
            "n_rows": int(self._n_rows),
            "n_index_keys": int(len(self.raw_index)),
            "nonempty_context_pools": int(nonempty),
            "mappable_keys_to_probs": int(ok_map),
            "probs_len": int(probs_len),
            "first_five_keys": list(list(self.raw_index.keys())[:5]),
        }
        with open(out_json, "w") as f:
            json.dump(stats, f, indent=2)

    def sample(
        self,
        signature_probs: np.ndarray,
        n: int,
        rng: np.random.Generator | None = None
    ) -> Tuple[pd.DataFrame, np.ndarray]:
        """Sample n mutations. We:
           - compute per-key probabilities only over keys that can be mapped,
           - renormalize to 1.0 over those keys,
           - draw multinomial counts over that subset,
           - sample rows from each pool (with replacement if needed)."""
        rng = rng or np.random.default_rng()
        probs_len = len(signature_probs)

        # Build key order and pvec over mappable keys that have non-empty pools
        key_order: List[Any] = []
        idxs: List[int] = []
        pvec: List[float] = []
        for k, pool in self.raw_index.items():
            if pool is None or len(pool) == 0:
                continue
            pos = self._key_to_prob_index(k, probs_len)
            if pos is None:
                continue
            key_order.append(k)
            idxs.append(pos)
            pvec.append(float(signature_probs[pos]))

        if not key_order:
            # nothing mappable / non-empty
            return (
                self.arrow_table.slice(0, 0).to_pandas(types_mapper=pd.ArrowDtype).copy(),
                np.array([], dtype=int),
            )

        pvec = np.asarray(pvec, dtype="float64")
        s = pvec.sum()
        if s <= 0:
            return (
                self.arrow_table.slice(0, 0).to_pandas(types_mapper=pd.ArrowDtype).copy(),
                np.array([], dtype=int),
            )
        pvec /= s

        draws = rng.multinomial(n, pvec)

        picked_rows: list[int] = []
        picked_ctx_positions: list[int] = []  # positions into probs (not the raw key)
        for k, pos, kdraw in zip(key_order, idxs, draws):
            if kdraw <= 0:
                continue
            pool = self.raw_index[k]
            replace = (kdraw > len(pool))
            rows = rng.choice(pool, size=kdraw, replace=replace)
            picked_rows.extend(rows.tolist())
            picked_ctx_positions.extend([pos] * len(rows))

        if not picked_rows:
            return (
                self.arrow_table.slice(0, 0).to_pandas(types_mapper=pd.ArrowDtype).copy(),
                np.array([], dtype=int),
            )

        df = self.arrow_table.take(np.asarray(picked_rows, dtype=np.int64)).to_pandas(
            types_mapper=pd.ArrowDtype
        ).copy()
        return df, np.asarray(picked_ctx_positions, dtype=int)
