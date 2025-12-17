#!/usr/bin/env python3
from __future__ import annotations

import argparse, json
from pathlib import Path
import pandas as pd
import numpy as np
import pysam

# ---------------- thresholds (NUMERIC ONLY) ----------------
SIFT4G_DAMAGING_LT = 0.05                 # <0.05 damaging, else benign

PP2_BENIGN_LT      = 0.15                 # PolyPhen-2 HDIV
PP2_POSS_UPPER     = 0.85                 # [0.15, 0.85) possibly damaging, >=0.85 damaging

CADD_BENIGN_LT     = 2.0                  # PHRED
CADD_UNCERT_LT     = 10.0                 # [2,10)
CADD_LIKELY_LT     = 20.0                 # [10,20), >=20 damaging

# ---------------- DBS78 list ----------------
DBS78 = [
    "AC>CA","AC>CG","AC>CT","AC>GA","AC>GG","AC>GT","AC>TA","AC>TG","AC>TT",
    "AT>CA","AT>CC","AT>CG","AT>GA","AT>GC","AT>TA",
    "CC>AA","CC>AG","CC>AT","CC>GA","CC>GG","CC>GT","CC>TA","CC>TG","CC>TT",
    "CG>AT","CG>GC","CG>GT","CG>TA","CG>TC","CG>TT",
    "CT>AA","CT>AC","CT>AG","CT>GA","CT>GC","CT>GG","CT>TA","CT>TC","CT>TG",
    "GC>AA","GC>AG","GC>AT","GC>CA","GC>CG","GC>TA",
    "TA>AT","TA>CG","TA>CT","TA>GC","TA>GG","TA>GT",
    "TC>AA","TC>AG","TC>AT","TC>CA","TC>CG","TC>CT","TC>GA","TC>GG","TC>GT",
    "TG>AA","TG>AC","TG>AT","TG>CA","TG>CC","TG>CT","TG>GA","TG>GC","TG>GT",
    "TT>AA","TT>AC","TT>AG","TT>CA","TT>CC","TT>CG","TT>GA","TT>GC","TT>GG",
]
ID2CTX = {i: s for i, s in enumerate(DBS78)}

DNA_COMP = str.maketrans("ACGT","TGCA")
def revcomp(s: str) -> str:
    return s.translate(DNA_COMP)[::-1]

def split_ctx(s: str) -> tuple[str,str]:
    a, b = s.split(">")
    return a, b

# ---------------- translation ----------------
CODON2AA = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L","CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M","GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S","CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T","GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*","CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K","GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W","CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R","GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}
def aa_of(codon: str) -> str|None:
    return CODON2AA.get(codon.upper().replace("U","T"))

def classify(ref_aa: str, alt_aa: str|None) -> str:
    if alt_aa is None:
        return "synonymous"
    if alt_aa == "*":
        return "stop_gained"
    return "missense" if alt_aa != ref_aa else "synonymous"

def worse(a: str, b: str) -> str:
    rank = {"stop_gained":0, "missense":1, "synonymous":2}
    return a if rank[a] <= rank[b] else b

# ---------------- score categories ----------------
def sift_cat(score: float|None) -> str:
    if score is None or np.isnan(score): return "missing"
    return "damaging" if score < SIFT4G_DAMAGING_LT else "benign"

def pp2_cat(score: float|None) -> str:
    if score is None or np.isnan(score): return "missing"
    if score < PP2_BENIGN_LT:  return "benign"
    if score < PP2_POSS_UPPER: return "possibly_damaging"
    return "damaging"

def cadd_cat(score: float|None) -> str:
    if score is None or np.isnan(score): return "missing"
    if score < CADD_BENIGN_LT:  return "benign"
    if score < CADD_UNCERT_LT:  return "uncertain"
    if score < CADD_LIKELY_LT:  return "likely_damaging"
    return "damaging"

# ---------------- dbNSFP AA-index loading ----------------
def _pick_col(cols: list[str], candidates: list[str]) -> str|None:
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None

def load_dbnsfp_aa_index(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)

    # required
    c_aaref = _pick_col(df.columns.tolist(), ["aaref"])
    c_aaalt = _pick_col(df.columns.tolist(), ["aaalt"])
    if c_aaref is None or c_aaalt is None:
        raise SystemExit(f"❌ dbNSFP AA index missing aaref/aaalt columns: {path}")

    # optional
    c_aapos = _pick_col(df.columns.tolist(), ["aapos"])

    # scores (try common names)
    c_sift = _pick_col(df.columns.tolist(), ["sift4g_score", "SIFT4G_score"])
    c_pp2  = _pick_col(df.columns.tolist(), ["pp2_hdiv_score", "Polyphen2_HDIV_score"])
    c_cadd = _pick_col(df.columns.tolist(), ["cadd_phred", "CADD_phred", "CADD_PHRED"])

    keep = [c_aaref, c_aaalt]
    if c_aapos: keep.append(c_aapos)
    for c in (c_sift, c_pp2, c_cadd):
        if c: keep.append(c)

    df = df[keep].copy()
    df = df.rename(columns={
        c_aaref: "aaref",
        c_aaalt: "aaalt",
        (c_aapos or "aapos"): "aapos",
        (c_sift or "sift4g_score"): "sift4g_score",
        (c_pp2 or "pp2_hdiv_score"): "pp2_hdiv_score",
        (c_cadd or "cadd_phred"): "cadd_phred",
    })

    # ensure numeric
    for sc in ["sift4g_score","pp2_hdiv_score","cadd_phred"]:
        if sc in df.columns:
            df[sc] = pd.to_numeric(df[sc], errors="coerce")

    return df

def build_score_lookup(db: pd.DataFrame, use_aapos: bool):
    """
    Returns two dicts:
      exact[(aaref, aaalt, aapos)] -> (sift, pp2, cadd)
      fallback[(aaref, aaalt)]     -> (sift, pp2, cadd)
    We aggregate per key using:
      - SIFT4G: min (lower worse)
      - PP2: max (higher worse)
      - CADD: max (higher worse)
    """
    def agg_scores(g: pd.DataFrame):
        sift = g["sift4g_score"].min(skipna=True) if "sift4g_score" in g else np.nan
        pp2  = g["pp2_hdiv_score"].max(skipna=True) if "pp2_hdiv_score" in g else np.nan
        cadd = g["cadd_phred"].max(skipna=True) if "cadd_phred" in g else np.nan
        return (None if np.isnan(sift) else float(sift),
                None if np.isnan(pp2)  else float(pp2),
                None if np.isnan(cadd) else float(cadd))

    fallback = {}
    for (aaref, aaalt), g in db.groupby(["aaref","aaalt"], sort=False):
        fallback[(aaref, aaalt)] = agg_scores(g)

    exact = {}
    if use_aapos and "aapos" in db.columns:
        for (aaref, aaalt, aapos), g in db.groupby(["aaref","aaalt","aapos"], sort=False):
            exact[(aaref, aaalt, int(aapos))] = agg_scores(g)

    return exact, fallback

# ---------------- AA consequence / AA change for one DBS row ----------------
def compute_primary_secondary(row, fasta: pysam.FastaFile):
    chrom = str(row["chr"])
    pos   = int(row["pos"])
    strand = int(row["strand"])
    ref_codon = str(row["ref_codon"]).upper()
    ref_aa = str(row["ref_aa"]).upper()
    ci = int(row["codon_index"])
    ctx = ID2CTX[int(row["context_id"])]
    _, alt2 = split_ctx(ctx)  # alt2 is 2bp

    # within-codon
    if ci in (0,1):
        alt = list(ref_codon)
        alt[ci]   = alt2[0]
        alt[ci+1] = alt2[1]
        alt_codon = "".join(alt)
        alt_aa = aa_of(alt_codon)
        lbl = classify(ref_aa, alt_aa)
        # primary AA change
        aaref_p, aaalt_p = ref_aa, (alt_aa if alt_aa is not None else ref_aa)
        return {
            "primary_label": lbl,
            "secondary_label": None,
            "primary_side": "single",
            "aaref_primary": aaref_p,
            "aaalt_primary": aaalt_p,
            "aaref_secondary": None,
            "aaalt_secondary": None,
        }

    # boundary: need the adjacent codon (transcript strand-aware)
    if strand == 1:
        next_start = pos + 1
        next_trip  = fasta.fetch(chrom, next_start-1, next_start-1+3).upper()
    else:
        next_start = pos - 1
        seq = fasta.fetch(chrom, next_start-3, next_start).upper()
        next_trip = revcomp(seq)

    ref_right_aa = aa_of(next_trip) or "X"

    # mutate last base of left codon + first base of right codon
    alt_left = list(ref_codon)
    alt_left[2] = alt2[0]
    alt_left = "".join(alt_left)
    alt_left_aa = aa_of(alt_left)

    alt_right = list(next_trip)
    alt_right[0] = alt2[1]
    alt_right = "".join(alt_right)
    alt_right_aa = aa_of(alt_right)

    lbl_left  = classify(ref_aa, alt_left_aa)
    lbl_right = classify(ref_right_aa, alt_right_aa)

    primary = worse(lbl_left, lbl_right)
    primary_side = "left" if primary == lbl_left else "right"
    secondary = lbl_right if primary_side == "left" else lbl_left

    if primary_side == "left":
        aaref_p, aaalt_p = ref_aa, (alt_left_aa if alt_left_aa is not None else ref_aa)
        aaref_s, aaalt_s = ref_right_aa, (alt_right_aa if alt_right_aa is not None else ref_right_aa)
    else:
        aaref_p, aaalt_p = ref_right_aa, (alt_right_aa if alt_right_aa is not None else ref_right_aa)
        aaref_s, aaalt_s = ref_aa, (alt_left_aa if alt_left_aa is not None else ref_aa)

    return {
        "primary_label": primary,
        "secondary_label": secondary,
        "primary_side": primary_side,
        "aaref_primary": aaref_p,
        "aaalt_primary": aaalt_p,
        "aaref_secondary": aaref_s,
        "aaalt_secondary": aaalt_s,
    }

def get_aapos_for_primary(row, primary_side: str):
    # If your parquet has aapos_primary/aapos_secondary, use them.
    # Otherwise return None and we’ll match by (aaref, aaalt) only.
    for cand in ["aapos_primary", "aa_pos_primary", "aapos"]:
        if cand in row and pd.notna(row[cand]):
            return int(row[cand])
    if primary_side == "right":
        for cand in ["aapos_secondary", "aa_pos_secondary"]:
            if cand in row and pd.notna(row[cand]):
                return int(row[cand])
    return None

# ---------------- main annotation ----------------
def main():
    ap = argparse.ArgumentParser(description="Annotate DBS sims with dbNSFP scores by matching AA changes (aaref/aaalt [+aapos]).")
    ap.add_argument("--parquet-in", required=True, help="rep_*.sim.with_aapos.parquet (or sim parquet without aapos, still works)")
    ap.add_argument("--dbnsfp-aa-index", required=True, help="ref/dbnsfp_aa_index.parquet (must include aaref, aaalt)")
    ap.add_argument("--fasta", default="ref/GRCh38.primary_assembly.genome.fa", help="Reference FASTA (indexed .fai) for boundary codon retrieval")
    ap.add_argument("--out-prefix", required=True, help="Output prefix (will write <prefix>_<rep>.json and optional TSV)")
    ap.add_argument("--write-tsv", action="store_true", help="Also write per-variant TSV with primary scores")
    args = ap.parse_args()

    pq_path = Path(args.parquet_in)
    if not pq_path.exists():
        raise SystemExit(f"❌ parquet not found: {pq_path}")

    fasta = pysam.FastaFile(args.fasta)

    # Load sim parquet (only columns that exist)
    # IMPORTANT: do NOT request aaref_primary, etc. Those are computed.
    sim_cols = ["chr","pos","strand","ref_codon","ref_aa","codon_index","context_id"]
    df = pd.read_parquet(pq_path, columns=[c for c in sim_cols if c in pd.read_parquet(pq_path, engine="pyarrow").columns])
    # If the above “columns exists” check is annoying, just read and slice:
    # df = pd.read_parquet(pq_path)[sim_cols]

    # Safer:
    all_df = pd.read_parquet(pq_path)
    missing = [c for c in sim_cols if c not in all_df.columns]
    if missing:
        raise SystemExit(f"❌ sim parquet missing required columns: {missing}\nFound: {list(all_df.columns)}")
    df = all_df[sim_cols + [c for c in all_df.columns if c.startswith("aapos_") or c.startswith("aa_pos_")]].copy()

    # Load dbNSFP AA index once and build lookup dicts
    db = load_dbnsfp_aa_index(Path(args.dbnsfp_aa_index))
    use_aapos = any(c in df.columns for c in ["aapos_primary","aa_pos_primary","aapos","aapos_secondary","aa_pos_secondary"])
    exact, fallback = build_score_lookup(db, use_aapos=use_aapos)

    # Precompute AA-change keys for all variants
    rows = []
    for _, r in df.iterrows():
        info = compute_primary_secondary(r, fasta)
        aapos_p = get_aapos_for_primary(r, info["primary_side"]) if use_aapos else None
        rows.append({
            "chr": r["chr"],
            "pos": int(r["pos"]),
            "codon_index": int(r["codon_index"]),
            **info,
            "aapos_primary": aapos_p,
        })
    ann = pd.DataFrame(rows)

    # Unique lookups only (FAST)
    uniq = ann[["aaref_primary","aaalt_primary","aapos_primary"]].drop_duplicates()

    score_cache = {}
    for _, u in uniq.iterrows():
        k2 = (u["aaref_primary"], u["aaalt_primary"])
        k3 = None
        if use_aapos and pd.notna(u["aapos_primary"]):
            k3 = (u["aaref_primary"], u["aaalt_primary"], int(u["aapos_primary"]))

        if k3 is not None and k3 in exact:
            score_cache[(k2, int(u["aapos_primary"]))] = exact[k3]
        else:
            score_cache[(k2, None if not use_aapos else (int(u["aapos_primary"]) if pd.notna(u["aapos_primary"]) else None))] = fallback.get(k2, (None, None, None))

    # Attach scores
    sift, pp2, cadd = [], [], []
    for _, r in ann.iterrows():
        k2 = (r["aaref_primary"], r["aaalt_primary"])
        k_aapos = int(r["aapos_primary"]) if (use_aapos and pd.notna(r["aapos_primary"])) else None
        v = score_cache.get((k2, k_aapos))
        if v is None:
            v = fallback.get(k2, (None, None, None))
        sift.append(v[0]); pp2.append(v[1]); cadd.append(v[2])

    ann["sift4g_score"] = sift
    ann["pp2_hdiv_score"] = pp2
    ann["cadd_phred"] = cadd

    # Summaries
    n = len(ann)
    sift_counts = {"damaging":0, "benign":0, "missing":0}
    pp2_counts  = {"benign":0, "possibly_damaging":0, "damaging":0, "missing":0}
    cadd_counts = {"benign":0, "uncertain":0, "likely_damaging":0, "damaging":0, "missing":0}

    for s in ann["sift4g_score"].tolist():
        sift_counts[sift_cat(s)] += 1
    for s in ann["pp2_hdiv_score"].tolist():
        pp2_counts[pp2_cat(s)] += 1
    for s in ann["cadd_phred"].tolist():
        cadd_counts[cadd_cat(s)] += 1

    rep_id = pq_path.stem
    out_prefix = Path(args.out_prefix)

    summary = {
        "replicate": rep_id,
        "n_variants": int(n),
        "SIFT4G": {
            "counts": sift_counts,
            "fractions": {k: v/max(n,1) for k,v in sift_counts.items()},
            "thresholds": {"damaging": f"< {SIFT4G_DAMAGING_LT}", "benign": f">= {SIFT4G_DAMAGING_LT}"}
        },
        "PolyPhen2_HDIV": {
            "counts": pp2_counts,
            "fractions": {k: v/max(n,1) for k,v in pp2_counts.items()},
            "thresholds": {
                "benign": f"< {PP2_BENIGN_LT}",
                "possibly_damaging": f"[{PP2_BENIGN_LT}, {PP2_POSS_UPPER})",
                "damaging": f">= {PP2_POSS_UPPER}"
            }
        },
        "CADD_PHRED": {
            "counts": cadd_counts,
            "fractions": {k: v/max(n,1) for k,v in cadd_counts.items()},
            "thresholds": {
                "benign": f"< {CADD_BENIGN_LT}",
                "uncertain": f"[{CADD_BENIGN_LT}, {CADD_UNCERT_LT})",
                "likely_damaging": f"[{CADD_UNCERT_LT}, {CADD_LIKELY_LT})",
                "damaging": f">= {CADD_LIKELY_LT}"
            }
        }
    }

    out_json = out_prefix.parent / f"{out_prefix.name}_{rep_id}.json"
    out_json.write_text(json.dumps(summary, indent=2))

    if args.write_tsv:
        out_tsv = out_prefix.parent / f"{out_prefix.name}_primary_scores_{rep_id}.tsv"
        ann_out = ann[[
            "chr","pos","codon_index",
            "primary_label","secondary_label","primary_side",
            "aaref_primary","aaalt_primary","aapos_primary",
            "sift4g_score","pp2_hdiv_score","cadd_phred"
        ]]
        ann_out.to_csv(out_tsv, sep="\t", index=False)

    print(f"✅ wrote {out_json}")
    if args.write_tsv:
        print(f"✅ wrote {out_tsv}")

if __name__ == "__main__":
    main()
