#!/usr/bin/env python3
from __future__ import annotations

import argparse, json
from pathlib import Path
from collections import Counter

import numpy as np
import pandas as pd
import pysam

# ---------------- thresholds (NUMERIC ONLY) ----------------
SIFT4G_DAMAGING_LT = 0.05

PP2_BENIGN_LT      = 0.15
PP2_POSS_UPPER     = 0.85

CADD_BENIGN_LT     = 2.0
CADD_UNCERT_LT     = 10.0
CADD_LIKELY_LT     = 20.0

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

DNA_COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s: str) -> str:
    return s.translate(DNA_COMP)[::-1]

def split_ctx(s: str) -> tuple[str, str]:
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

def aa_of(codon: str) -> str | None:
    codon = str(codon).upper().replace("U", "T")
    return CODON2AA.get(codon)

def classify(ref_aa: str, alt_aa: str | None) -> str:
    if alt_aa is None:
        return "synonymous"
    if alt_aa == "*":
        return "stop_gained"
    return "missense" if alt_aa != ref_aa else "synonymous"

def worse(a: str, b: str) -> str:
    rank = {"stop_gained": 0, "missense": 1, "synonymous": 2}
    return a if rank[a] <= rank[b] else b

# ---------------- score categories ----------------
def sift_cat(score: float | None) -> str:
    if score is None or (isinstance(score, float) and np.isnan(score)): return "missing"
    return "damaging" if score < SIFT4G_DAMAGING_LT else "benign"

def pp2_cat(score: float | None) -> str:
    if score is None or (isinstance(score, float) and np.isnan(score)): return "missing"
    if score < PP2_BENIGN_LT:  return "benign"
    if score < PP2_POSS_UPPER: return "possibly_damaging"
    return "damaging"

def cadd_cat(score: float | None) -> str:
    if score is None or (isinstance(score, float) and np.isnan(score)): return "missing"
    if score < CADD_BENIGN_LT:  return "benign"
    if score < CADD_UNCERT_LT:  return "uncertain"
    if score < CADD_LIKELY_LT:  return "likely_damaging"
    return "damaging"

# ---------------- dbNSFP AA-index loading ----------------
def _pick_col(cols: list[str], candidates: list[str]) -> str | None:
    low = {c.lower(): c for c in cols}
    for cand in candidates:
        if cand.lower() in low:
            return low[cand.lower()]
    return None

def load_dbnsfp_aa_index(path: Path) -> pd.DataFrame:
    df = pd.read_parquet(path)
    cols = df.columns.tolist()

    c_aaref = _pick_col(cols, ["aaref"])
    c_aaalt = _pick_col(cols, ["aaalt"])
    if c_aaref is None or c_aaalt is None:
        raise SystemExit(f"❌ dbNSFP AA index missing aaref/aaalt columns: {path}")

    c_aapos = _pick_col(cols, ["aapos"])
    c_sift  = _pick_col(cols, ["sift4g_score", "SIFT4G_score"])
    c_pp2   = _pick_col(cols, ["pp2_hdiv_score", "Polyphen2_HDIV_score"])
    c_cadd  = _pick_col(cols, ["cadd_phred", "CADD_phred", "CADD_PHRED"])

    keep = [c_aaref, c_aaalt]
    if c_aapos: keep.append(c_aapos)
    if c_sift:  keep.append(c_sift)
    if c_pp2:   keep.append(c_pp2)
    if c_cadd:  keep.append(c_cadd)

    df = df[keep].copy()
    ren = {c_aaref: "aaref", c_aaalt: "aaalt"}
    if c_aapos: ren[c_aapos] = "aapos"
    if c_sift:  ren[c_sift]  = "sift4g_score"
    if c_pp2:   ren[c_pp2]   = "pp2_hdiv_score"
    if c_cadd:  ren[c_cadd]  = "cadd_phred"
    df = df.rename(columns=ren)

    for sc in ["sift4g_score", "pp2_hdiv_score", "cadd_phred"]:
        if sc in df.columns:
            df[sc] = pd.to_numeric(df[sc], errors="coerce")

    if "aapos" in df.columns:
        df["aapos"] = pd.to_numeric(df["aapos"], errors="coerce")

    # normalize AA strings
    df["aaref"] = df["aaref"].astype(str).str.strip()
    df["aaalt"] = df["aaalt"].astype(str).str.strip()

    return df

def build_score_lookup(db: pd.DataFrame, use_aapos: bool):
    """
    exact[(aaref, aaalt, aapos)] -> (sift_min, pp2_max, cadd_max)
    fallback[(aaref, aaalt)]     -> (sift_min, pp2_max, cadd_max)
    """
    def agg_scores(g: pd.DataFrame):
        sift = g["sift4g_score"].min(skipna=True) if "sift4g_score" in g.columns else np.nan
        pp2  = g["pp2_hdiv_score"].max(skipna=True) if "pp2_hdiv_score" in g.columns else np.nan
        cadd = g["cadd_phred"].max(skipna=True) if "cadd_phred" in g.columns else np.nan
        return (
            None if (not np.isfinite(sift)) else float(sift),
            None if (not np.isfinite(pp2))  else float(pp2),
            None if (not np.isfinite(cadd)) else float(cadd),
        )

    fallback = {}
    for (aaref, aaalt), g in db.groupby(["aaref", "aaalt"], sort=False):
        fallback[(aaref, aaalt)] = agg_scores(g)

    exact = {}
    if use_aapos and "aapos" in db.columns:
        db2 = db.dropna(subset=["aapos"]).copy()
        db2["aapos"] = db2["aapos"].astype(int)
        for (aaref, aaalt, aapos), g in db2.groupby(["aaref", "aaalt", "aapos"], sort=False):
            exact[(aaref, aaalt, int(aapos))] = agg_scores(g)

    return exact, fallback

# ---------------- AA consequence / AA change for one DBS row ----------------
def compute_primary_secondary(row, fasta: pysam.FastaFile):
    chrom     = str(row["chr"])
    pos       = int(row["pos"])       # 1-based
    strand    = int(row["strand"])    # +1 / -1
    ref_codon = str(row["ref_codon"]).upper()
    ref_aa    = str(row["ref_aa"]).upper()
    ctx       = ID2CTX[int(row["context_id"])]

    # Use codon_pos (0/1/2) if present, otherwise codon_index (old meaning).
    if "codon_pos" in row and pd.notna(row["codon_pos"]):
        ci = int(row["codon_pos"])
    else:
        ci = int(row["codon_index"])

    _, alt2 = split_ctx(ctx)  # alt2 always written forward in your context list

    # IMPORTANT: if the transcript strand is -1, the actual ALT applied on genome is revcomp(alt2)
    # because the simulation "strand" indicates coding strand direction.
    if strand == -1:
        alt2 = revcomp(alt2)

    # within-codon (codon_pos 0 or 1)
    if ci in (0, 1):
        alt = list(ref_codon)
        alt[ci]   = alt2[0]
        alt[ci+1] = alt2[1]
        alt_codon = "".join(alt)
        alt_aa = aa_of(alt_codon)
        lbl = classify(ref_aa, alt_aa)
        return {
            "primary_label": lbl,
            "secondary_label": None,
            "primary_side": "single",
            "aaref_primary": ref_aa,
            "aaalt_primary": (alt_aa if alt_aa is not None else ref_aa),
            "aaref_secondary": None,
            "aaalt_secondary": None,
        }

    # boundary (codon_pos 2): left codon gets base0, right codon gets base1
    if strand == 1:
        next_start = pos + 1  # next base after the DBS start
        next_trip  = fasta.fetch(chrom, next_start-1, next_start-1+3).upper()
    else:
        next_start = pos - 1
        seq = fasta.fetch(chrom, next_start-3, next_start).upper()
        next_trip = revcomp(seq)

    ref_right_aa = aa_of(next_trip) or "X"

    alt_left = list(ref_codon)
    alt_left[2] = alt2[0]
    alt_left = "".join(alt_left)
    alt_left_aa = aa_of(alt_left)

    if len(next_trip) == 3 and "N" not in next_trip:
        alt_right = list(next_trip)
        alt_right[0] = alt2[1]
        alt_right = "".join(alt_right)
        alt_right_aa = aa_of(alt_right)
    else:
        alt_right_aa = None

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

def detect_sim_aapos_columns(df: pd.DataFrame) -> list[str]:
    # accept a few possibilities if you created these in older pipelines
    candidates = [
        "aapos_primary", "aa_pos_primary",
        "aapos", "aa_pos",
        "aapos_secondary", "aa_pos_secondary",
    ]
    return [c for c in candidates if c in df.columns]

def get_aapos_for_primary(row: pd.Series, primary_side: str) -> int | None:
    # best-effort; if you do NOT have aapos, we will just return None
    for cand in ["aapos_primary", "aa_pos_primary", "aapos", "aa_pos"]:
        if cand in row.index and pd.notna(row[cand]):
            return int(row[cand])
    if primary_side == "right":
        for cand in ["aapos_secondary", "aa_pos_secondary"]:
            if cand in row.index and pd.notna(row[cand]):
                return int(row[cand])
    return None

# ---------------- main ----------------
def main():
    ap = argparse.ArgumentParser(
        description="Annotate DBS sims with dbNSFP scores by matching AA changes (aaref/aaalt [+aapos if present])."
    )
    ap.add_argument("--parquet-in", required=True, help="rep_*.sim.parquet (new) OR rep_*.sim.with_aapos.parquet (old)")
    ap.add_argument("--dbnsfp-aa-index", required=True, help="ref/dbnsfp_aa_index.parquet (must include aaref, aaalt)")
    ap.add_argument("--fasta", default="ref/GRCh38.primary_assembly.genome.fa", help="Reference FASTA (indexed .fai)")
    ap.add_argument("--out-prefix", required=True, help="Output prefix (writes <prefix>_<rep>.json and optional TSV)")
    ap.add_argument("--write-tsv", action="store_true", help="Also write per-variant TSV with primary scores")
    args = ap.parse_args()

    pq_path = Path(args.parquet_in)
    if not pq_path.exists():
        raise SystemExit(f"❌ parquet not found: {pq_path}")

    fasta = pysam.FastaFile(args.fasta)

    sim_required = ["chr","pos","strand","ref_codon","ref_aa","codon_index","context_id"]
    all_df = pd.read_parquet(pq_path)
    missing = [c for c in sim_required if c not in all_df.columns]
    if missing:
        raise SystemExit(f"❌ sim parquet missing required columns: {missing}\nFound: {list(all_df.columns)}")

    # keep required + optional codon_pos + any aapos-like cols
    keep = sim_required.copy()
    if "codon_pos" in all_df.columns:
        keep.append("codon_pos")
    keep += detect_sim_aapos_columns(all_df)
    df = all_df[keep].copy()

    # Load dbNSFP AA-index once and build lookup dicts
    db = load_dbnsfp_aa_index(Path(args.dbnsfp_aa_index))
    have_aapos = len(detect_sim_aapos_columns(df)) > 0
    exact, fallback = build_score_lookup(db, use_aapos=have_aapos)

    # Compute AA-change + primary/secondary for each row
    out_rows = []
    for _, r in df.iterrows():
        info = compute_primary_secondary(r, fasta)
        aapos_p = get_aapos_for_primary(r, info["primary_side"]) if have_aapos else None
        out_rows.append({
            "chr": r["chr"],
            "pos": int(r["pos"]),
            "codon_index": int(r["codon_index"]),
            "codon_pos": (int(r["codon_pos"]) if "codon_pos" in r.index and pd.notna(r["codon_pos"]) else None),
            **info,
            "aapos_primary": aapos_p,
        })
    ann = pd.DataFrame(out_rows)

    # Lookup scores (unique keys only)
    uniq = ann[["aaref_primary","aaalt_primary","aapos_primary"]].drop_duplicates()

    cache = {}
    for _, u in uniq.iterrows():
        aaref = str(u["aaref_primary"])
        aaalt = str(u["aaalt_primary"])
        k2 = (aaref, aaalt)
        if have_aapos and pd.notna(u["aapos_primary"]):
            k3 = (aaref, aaalt, int(u["aapos_primary"]))
            if k3 in exact:
                cache[(k2, int(u["aapos_primary"]))] = exact[k3]
                continue
        cache[(k2, None)] = fallback.get(k2, (None, None, None))

    sift, pp2, cadd = [], [], []
    for _, r in ann.iterrows():
        k2 = (str(r["aaref_primary"]), str(r["aaalt_primary"]))
        key = (k2, int(r["aapos_primary"])) if (have_aapos and pd.notna(r["aapos_primary"])) else (k2, None)
        v = cache.get(key)
        if v is None:
            v = fallback.get(k2, (None, None, None))
        sift.append(v[0]); pp2.append(v[1]); cadd.append(v[2])

    ann["sift4g_score"] = sift
    ann["pp2_hdiv_score"] = pp2
    ann["cadd_phred"] = cadd

    # Summaries
    n = len(ann)
    sift_counts = Counter(sift_cat(x) for x in ann["sift4g_score"].tolist())
    pp2_counts  = Counter(pp2_cat(x) for x in ann["pp2_hdiv_score"].tolist())
    cadd_counts = Counter(cadd_cat(x) for x in ann["cadd_phred"].tolist())

    rep_id = pq_path.stem
    out_prefix = Path(args.out_prefix)
    out_dir = out_prefix.parent
    out_dir.mkdir(parents=True, exist_ok=True)

    summary = {
        "replicate": rep_id,
        "n_variants": int(n),
        "SIFT4G": {
            "counts": dict(sift_counts),
            "fractions": {k: float(v)/max(n,1) for k,v in sift_counts.items()},
            "thresholds": {"damaging": f"< {SIFT4G_DAMAGING_LT}", "benign": f">= {SIFT4G_DAMAGING_LT}"}
        },
        "PolyPhen2_HDIV": {
            "counts": dict(pp2_counts),
            "fractions": {k: float(v)/max(n,1) for k,v in pp2_counts.items()},
            "thresholds": {
                "benign": f"< {PP2_BENIGN_LT}",
                "possibly_damaging": f"[{PP2_BENIGN_LT}, {PP2_POSS_UPPER})",
                "damaging": f">= {PP2_POSS_UPPER}"
            }
        },
        "CADD_PHRED": {
            "counts": dict(cadd_counts),
            "fractions": {k: float(v)/max(n,1) for k,v in cadd_counts.items()},
            "thresholds": {
                "benign": f"< {CADD_BENIGN_LT}",
                "uncertain": f"[{CADD_BENIGN_LT}, {CADD_UNCERT_LT})",
                "likely_damaging": f"[{CADD_UNCERT_LT}, {CADD_LIKELY_LT})",
                "damaging": f">= {CADD_LIKELY_LT}"
            }
        },
        "matching_mode": "aaref+aaalt+aapos" if have_aapos else "aaref+aaalt",
        "inputs": {
            "parquet": str(pq_path),
            "dbnsfp_aa_index": str(Path(args.dbnsfp_aa_index)),
            "fasta": str(Path(args.fasta)),
        }
    }

    out_json = out_dir / f"{out_prefix.name}_{rep_id}.json"
    out_json.write_text(json.dumps(summary, indent=2))

    if args.write_tsv:
        out_tsv = out_dir / f"{out_prefix.name}_primary_scores_{rep_id}.tsv"
        ann_out = ann[[
            "chr","pos","codon_index","codon_pos",
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
