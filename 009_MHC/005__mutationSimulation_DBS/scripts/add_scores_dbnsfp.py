#!/usr/bin/env python3
"""
Annotate DBS simulations at AA/codon level and summarize pathogenicity scores
using dbNSFP AA fields (aaref/aaalt[/aapos]) as a proxy.

Key idea:
- dbNSFP is queried by genomic SNVs (chr/pos/ref/alt),
- but we KEEP only rows whose (aaref, aaalt, and optionally aapos) match the
  simulated AA change for the codon side we are evaluating.

Outputs (per replicate):
- JSON summary with counts/fractions per score class
- Optional TSV with per-variant primary/secondary labels and matched scores

Requirements:
- pandas, pysam
- tabix available in PATH
- dbNSFP bgzip+tabix indexed
"""

import argparse
import json
import re
import subprocess
from dataclasses import dataclass
from pathlib import Path
from typing import Optional, Dict, Tuple, List

import pandas as pd
import pysam

# ---------------- Codon table (standard) ----------------
CODON2AA = {
    # Phenylalanine / Leucine
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    # Isoleucine / Methionine
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    # Valine
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    # Serine
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "AGT":"S","AGC":"S",
    # Proline
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    # Threonine
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    # Alanine
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    # Tyrosine / STOP
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    # Histidine / Glutamine
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    # Asparagine / Lysine
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    # Aspartate / Glutamate
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    # Cysteine / Tryptophan / STOP
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    # Arginine
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGA":"R","AGG":"R",
    # Glycine
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G",
}

BASES = ["A","C","G","T"]
COMP = str.maketrans("ACGTacgt", "TGCAtgca")

def revcomp(seq: str) -> str:
    return seq.translate(COMP)[::-1]

# ---------------- thresholds from your screenshot ----------------
def classify_sift4g(score: Optional[float]) -> str:
    if score is None:
        return "missing"
    return "damaging" if score < 0.05 else "benign"

def classify_pp2_hdiv(score: Optional[float]) -> str:
    if score is None:
        return "missing"
    if score < 0.15:
        return "benign"
    if score < 0.85:
        return "possibly_damaging"
    return "damaging"

def classify_cadd(phred: Optional[float]) -> str:
    if phred is None:
        return "missing"
    if phred < 2:
        return "benign"
    if phred < 10:
        return "uncertain"
    if phred < 20:
        return "likely_damaging"
    return "damaging"

# ---------------- dbNSFP helper ----------------
@dataclass
class DBNSFPCols:
    IDX_CHR: int = 1-1
    IDX_POS: int = 2-1
    IDX_REF: int = 3-1
    IDX_ALT: int = 4-1
    IDX_AAREF: int = 5-1
    IDX_AAALT: int = 6-1
    IDX_AAPOS: int = 12-1

    IDX_TRANSCRIPTID: int = 15-1
    IDX_VEP_CANONICAL: int = 26-1

    IDX_SIFT4G_SCORE: int = 50-1
    IDX_PP2_HDIV_SCORE: int = 53-1
    IDX_CADD_PHRED: int = 144-1

DB = DBNSFPCols()

def tabix_fetch(path: Path, chrom: str, pos: int) -> List[List[str]]:
    # try both chr styles
    chroms = [chrom, chrom.lstrip("chr"), f"chr{chrom.lstrip('chr')}"]
    seen = set()
    for c in chroms:
        if c in seen:
            continue
        seen.add(c)
        try:
            res = subprocess.run(
                ["tabix", str(path), f"{c}:{pos}-{pos}"],
                capture_output=True, check=True
            )
            out = res.stdout.decode("utf-8", errors="replace").strip()
            if not out:
                continue
            rows = [ln.split("\t") for ln in out.split("\n") if ln.strip()]
            if rows:
                return rows
        except subprocess.CalledProcessError:
            pass
    return []

def split_or_none(s: str) -> Optional[List[str]]:
    if not s or s == ".":
        return None
    return s.split(";")

def get_canonical_index(canon_flags: Optional[List[str]]) -> Optional[int]:
    if not canon_flags:
        return None
    try:
        return canon_flags.index("YES")
    except ValueError:
        return None

def safe_float(x: Optional[str]) -> Optional[float]:
    if x is None or x == "." or x == "":
        return None
    try:
        return float(x)
    except ValueError:
        return None

def pick_score_at_index(values: Optional[List[str]], idx: Optional[int]) -> Optional[float]:
    if values and idx is not None and 0 <= idx < len(values):
        return safe_float(values[idx])
    # fallback: first numeric
    if values:
        for v in values:
            f = safe_float(v)
            if f is not None:
                return f
    return None

def best_dbnsfp_match_for_aa_change(
    dbnsfp_path: Path,
    chrom: str,
    pos: int,
    ref: str,
    alt: str,
    aa_ref: str,
    aa_alt: str,
    aa_pos: Optional[int] = None,
) -> Dict[str, Optional[float]]:
    """
    Query dbNSFP at chr:pos with SNV ref/alt, then keep only rows matching the AA change.
    If aa_pos is provided, also require dbNSFP aapos match.
    Uses canonical transcript if available; otherwise falls back to first numeric.
    """
    rows = tabix_fetch(dbnsfp_path, chrom, pos)
    # filter by ref/alt first
    rows = [r for r in rows if len(r) > DB.IDX_CADD_PHRED and r[DB.IDX_REF] == ref and r[DB.IDX_ALT] == alt]
    if not rows:
        return {"sift4g_score": None, "pp2_hdiv_score": None, "cadd_phred": None}

    # filter by AA change (and position if available)
    def aa_ok(r: List[str]) -> bool:
        if r[DB.IDX_AAREF] != aa_ref or r[DB.IDX_AAALT] != aa_alt:
            return False
        if aa_pos is None:
            return True
        try:
            return int(r[DB.IDX_AAPOS]) == int(aa_pos)
        except:
            return False

    rows2 = [r for r in rows if aa_ok(r)]
    if not rows2:
        return {"sift4g_score": None, "pp2_hdiv_score": None, "cadd_phred": None}

    # choose canonical transcript if possible
    # dbNSFP columns are semicolon lists per transcript for many scores, but REF/ALT are per row
    # We'll pick using canonical index if present, else first numeric.
    r = rows2[0]
    canon_flags = split_or_none(r[DB.IDX_VEP_CANONICAL])
    idx = get_canonical_index(canon_flags)

    sift_scores = split_or_none(r[DB.IDX_SIFT4G_SCORE])
    pp2_scores  = split_or_none(r[DB.IDX_PP2_HDIV_SCORE])

    return {
        "sift4g_score": pick_score_at_index(sift_scores, idx),
        "pp2_hdiv_score": pick_score_at_index(pp2_scores, idx),
        "cadd_phred": safe_float(r[DB.IDX_CADD_PHRED]),
    }

# ---------------- DBS context mapping ----------------
def load_contexts_list(path: Path) -> List[str]:
    # file can be one per line OR tab/space separated list
    txt = path.read_text().strip()
    lines = [ln.strip() for ln in txt.splitlines() if ln.strip()]
    if len(lines) == 1 and ("\t" in lines[0] or " " in lines[0]):
        toks = re.split(r"[\t ]+", lines[0].strip())
        return [t for t in toks if t]
    return lines

def ctxid_to_ref2_alt2(ctx_list: List[str], ctx_id: int) -> Tuple[str, str]:
    # expects strings like "AC>CA"
    s = ctx_list[int(ctx_id)]
    left, right = s.split(">")
    return left, right

# ---------------- Codon / consequence logic ----------------
def translate_codon(c: str) -> str:
    c = c.upper()
    return CODON2AA.get(c, "X")  # X = unknown/invalid

def label_consequence(aa_ref: str, aa_alt: str) -> str:
    if aa_alt == aa_ref:
        return "synonymous"
    if aa_alt == "*" and aa_ref != "*":
        return "stop_gained"
    return "missense"

def worse_label(a: str, b: str) -> str:
    # STOP > missense > synonymous
    rank = {"synonymous": 0, "missense": 1, "stop_gained": 2}
    return a if rank.get(a, -1) >= rank.get(b, -1) else b

def mutate_single_codon(ref_codon: str, codon_index: int, ref2: str, alt2: str, strand: int) -> str:
    """
    Apply dinuc change to one codon (codon_index 0 or 1).
    We treat ref_codon as coding (transcript 5'->3').

    If strand == -1, ref2/alt2 are in genomic orientation; convert to transcript orientation by revcomp.
    """
    cod = ref_codon.upper()
    if len(cod) != 3:
        return cod

    if strand == -1:
        ref2_t = revcomp(ref2)
        alt2_t = revcomp(alt2)
    else:
        ref2_t = ref2
        alt2_t = alt2

    i = int(codon_index)
    # i=0 affects cod[0:2], i=1 affects cod[1:3]
    seg = cod[i:i+2]
    # if mismatch, still apply replacement (simulation might have already guaranteed, but keep robust)
    if seg != ref2_t:
        # do not block; replace anyway to stay consistent with "event happened"
        pass
    newcod = list(cod)
    newcod[i] = alt2_t[0]
    newcod[i+1] = alt2_t[1]
    return "".join(newcod)

def mutate_two_codons(ref_codon: str, next_ref_codon: str, ref2: str, alt2: str, strand: int) -> Tuple[str,str]:
    """
    codon_index==2: affects last base of this codon + first base of next codon.
    Returns (mutated_left_codon, mutated_right_codon)
    """
    left = ref_codon.upper()
    right = next_ref_codon.upper() if isinstance(next_ref_codon, str) else ""
    if len(left) != 3 or len(right) != 3:
        return left, right

    if strand == -1:
        ref2_t = revcomp(ref2)
        alt2_t = revcomp(alt2)
    else:
        ref2_t = ref2
        alt2_t = alt2

    # transcript orientation: left codon base3 + right codon base1
    # (index 2 in left, index 0 in right)
    # Again, tolerate mismatch
    new_left = list(left)
    new_right = list(right)
    new_left[2] = alt2_t[0]
    new_right[0] = alt2_t[1]
    return "".join(new_left), "".join(new_right)

# ---------------- Main per-parquet processing ----------------
def infer_next_codon(df: pd.DataFrame) -> pd.Series:
    """
    Best-effort: for codon_index==2 we need the 'next' codon.
    If your parquet already has something like next_ref_codon, we use it.
    Otherwise we approximate by shifting within transcript_id ordered by genomic pos.
    """
    for col in ["next_ref_codon", "ref_codon_next", "next_codon"]:
        if col in df.columns:
            return df[col].astype(str)

    # fallback: sort by transcript_id then by pos, shift ref_codon by -1
    if "transcript_id" in df.columns and "pos" in df.columns:
        tmp = df.sort_values(["transcript_id", "pos"]).copy()
        nxt = tmp.groupby("transcript_id")["ref_codon"].shift(-1)
        # map back
        nxt = nxt.reindex(df.index)
        return nxt.astype(str)

    return pd.Series([""]*len(df), index=df.index)

def get_aa_pos_column(df: pd.DataFrame) -> Optional[str]:
    for col in ["aa_pos", "aapos", "aa_position", "protein_pos", "ref_aapos"]:
        if col in df.columns:
            return col
    return None

def process_replicate(
    parquet_in: Path,
    contexts_list: List[str],
    dbnsfp: Path,
    out_json: Path,
    out_tsv: Optional[Path] = None,
):
    df = pd.read_parquet(parquet_in)
    need = ["chr","pos","strand","codon_index","ref_codon","ref_aa","context_id"]
    miss = [c for c in need if c not in df.columns]
    if miss:
        raise SystemExit(f"❌ Missing required columns in {parquet_in}: {miss}")

    aa_pos_col = get_aa_pos_column(df)
    if aa_pos_col is None:
        print("⚠️ No AA-position column found in parquet (aapos/aa_pos/etc). "
              "Will match dbNSFP proxies by aaref/aaalt only (less specific).")

    next_codon = infer_next_codon(df)

    # per-variant results (primary & secondary)
    rows = []
    for i, r in df.iterrows():
        chrom = str(r["chr"])
        pos = int(r["pos"])
        strand = int(r["strand"])
        ci = int(r["codon_index"])
        ref_codon = str(r["ref_codon"]).upper()
        ref_aa = str(r["ref_aa"]).upper()
        ctx_id = int(r["context_id"])

        ref2, alt2 = ctxid_to_ref2_alt2(contexts_list, ctx_id)

        # compute AA changes
        if ci in (0,1):
            alt_codon = mutate_single_codon(ref_codon, ci, ref2, alt2, strand)
            aa_alt = translate_codon(alt_codon)
            primary_label = label_consequence(ref_aa, aa_alt)

            secondary_label = None
            primary_side = "single"
            aa_ref_primary, aa_alt_primary = ref_aa, aa_alt
            aa_pos = int(r[aa_pos_col]) if aa_pos_col else None

            # enumerate proxy SNVs within codon: try each base position substitution
            # We must query dbNSFP by chr/pos/ref/alt (SNV) at the 2bp positions only is fine too,
            # but proxy should represent any SNV inside this codon that yields same AA change.
            # Minimal approach: query only the two mutated positions as SNVs using ref2->alt2 bases.
            # (This is usually enough to find AA-matching proxies; and keeps runtime sane.)
            # Position mapping: use DBS genomic pos and strand:
            # - the dinuc is at pos,pos+1 on the genome (1-based, increasing coords).
            # We make SNVs at those 2 positions.
            ref2g, alt2g = (ref2, alt2)  # genomic orientation already
            snv_candidates = [
                (pos,     ref2g[0], alt2g[0]),
                (pos + 1, ref2g[1], alt2g[1]),
            ]

            best = {"sift4g_score": None, "pp2_hdiv_score": None, "cadd_phred": None}
            for p_snv, refb, altb in snv_candidates:
                hit = best_dbnsfp_match_for_aa_change(
                    dbnsfp, chrom, p_snv, refb, altb, aa_ref_primary, aa_alt_primary, aa_pos=aa_pos
                )
                # keep the one with any non-missing score (preference: CADD then PP2 then SIFT)
                def score_key(d):
                    return (
                        d["cadd_phred"] is not None,
                        d["pp2_hdiv_score"] is not None,
                        d["sift4g_score"] is not None
                    )
                if score_key(hit) > score_key(best):
                    best = hit

            rows.append({
                "chr": chrom, "pos": pos, "codon_index": ci,
                "primary_label": primary_label,
                "secondary_label": secondary_label,
                "primary_side": primary_side,
                "sift4g_score": best["sift4g_score"],
                "pp2_hdiv_score": best["pp2_hdiv_score"],
                "cadd_phred": best["cadd_phred"],
            })

        elif ci == 2:
            # spans two codons
            ref_codon_R = str(next_codon.loc[i]).upper()
            altL, altR = mutate_two_codons(ref_codon, ref_codon_R, ref2, alt2, strand)

            aa_alt_L = translate_codon(altL)
            aa_alt_R = translate_codon(altR)

            label_L = label_consequence(ref_aa, aa_alt_L)
            # ref_aa for right codon may not be available; best-effort: translate ref_codon_R
            ref_aa_R = translate_codon(ref_codon_R)
            label_R = label_consequence(ref_aa_R, aa_alt_R)

            # primary = most deleterious label
            if worse_label(label_L, label_R) == label_L:
                primary_side = "left"
                primary_label = label_L
                secondary_side = "right"
                secondary_label = label_R
                aa_ref_primary, aa_alt_primary = ref_aa, aa_alt_L
                aa_ref_secondary, aa_alt_secondary = ref_aa_R, aa_alt_R
            else:
                primary_side = "right"
                primary_label = label_R
                secondary_side = "left"
                secondary_label = label_L
                aa_ref_primary, aa_alt_primary = ref_aa_R, aa_alt_R
                aa_ref_secondary, aa_alt_secondary = ref_aa, aa_alt_L

            aa_pos = int(r[aa_pos_col]) if aa_pos_col else None

            # proxy SNVs: two positions of the dinuc
            ref2g, alt2g = (ref2, alt2)
            snv_candidates = [
                (pos,     ref2g[0], alt2g[0]),
                (pos + 1, ref2g[1], alt2g[1]),
            ]

            # primary proxy
            bestP = {"sift4g_score": None, "pp2_hdiv_score": None, "cadd_phred": None}
            for p_snv, refb, altb in snv_candidates:
                hit = best_dbnsfp_match_for_aa_change(
                    dbnsfp, chrom, p_snv, refb, altb, aa_ref_primary, aa_alt_primary, aa_pos=aa_pos
                )
                def score_key(d):
                    return (d["cadd_phred"] is not None, d["pp2_hdiv_score"] is not None, d["sift4g_score"] is not None)
                if score_key(hit) > score_key(bestP):
                    bestP = hit

            # secondary proxy (optional, keep it for QC)
            bestS = {"sift4g_score": None, "pp2_hdiv_score": None, "cadd_phred": None}
            for p_snv, refb, altb in snv_candidates:
                hit = best_dbnsfp_match_for_aa_change(
                    dbnsfp, chrom, p_snv, refb, altb, aa_ref_secondary, aa_alt_secondary, aa_pos=aa_pos
                )
                def score_key(d):
                    return (d["cadd_phred"] is not None, d["pp2_hdiv_score"] is not None, d["sift4g_score"] is not None)
                if score_key(hit) > score_key(bestS):
                    bestS = hit

            rows.append({
                "chr": chrom, "pos": pos, "codon_index": ci,
                "primary_label": primary_label,
                "secondary_label": secondary_label,
                "primary_side": primary_side,
                "sift4g_score": bestP["sift4g_score"],
                "pp2_hdiv_score": bestP["pp2_hdiv_score"],
                "cadd_phred": bestP["cadd_phred"],
                "secondary_sift4g_score": bestS["sift4g_score"],
                "secondary_pp2_hdiv_score": bestS["pp2_hdiv_score"],
                "secondary_cadd_phred": bestS["cadd_phred"],
            })
        else:
            # unexpected codon_index
            continue

    out_df = pd.DataFrame(rows)
    n = len(out_df)

    # classify scores
    sift_cls = out_df["sift4g_score"].apply(classify_sift4g)
    pp2_cls  = out_df["pp2_hdiv_score"].apply(classify_pp2_hdiv)
    cadd_cls = out_df["cadd_phred"].apply(classify_cadd)

    def summarize_classes(series: pd.Series, order: List[str]) -> Dict:
        counts = {k: int((series == k).sum()) for k in order}
        counts["missing"] = int((series == "missing").sum()) if "missing" not in counts else counts["missing"]
        denom = max(n, 1)
        fracs = {k: round(counts[k]/denom, 6) for k in counts}
        return {"counts": counts, "fractions": fracs}

    summary = {
        "replicate": parquet_in.stem,
        "n_variants": int(n),
        "note": "Scores are dbNSFP proxy SNVs filtered by AA change (aaref/aaalt[/aapos]) to match simulated AA effects.",
        "thresholds": {
            "SIFT4G": {"damaging": "< 0.05", "benign": ">= 0.05"},
            "PolyPhen2_HDIV": {"benign": "< 0.15", "possibly_damaging": "0.15-0.85", "damaging": ">= 0.85"},
            "CADD_PHRED": {"benign": "< 2", "uncertain": "2-10", "likely_damaging": "10-20", "damaging": ">= 20"},
        },
        "SIFT4G": summarize_classes(sift_cls, ["damaging","benign","missing"]),
        "PolyPhen2_HDIV": summarize_classes(pp2_cls, ["benign","possibly_damaging","damaging","missing"]),
        "CADD_PHRED": summarize_classes(cadd_cls, ["benign","uncertain","likely_damaging","damaging","missing"]),
        "consequence_primary": {
            "counts": {
                "synonymous": int((out_df["primary_label"]=="synonymous").sum()),
                "missense": int((out_df["primary_label"]=="missense").sum()),
                "stop_gained": int((out_df["primary_label"]=="stop_gained").sum()),
            },
            "fractions": {
                "synonymous": round((out_df["primary_label"]=="synonymous").mean(), 6) if n else 0.0,
                "missense": round((out_df["primary_label"]=="missense").mean(), 6) if n else 0.0,
                "stop_gained": round((out_df["primary_label"]=="stop_gained").mean(), 6) if n else 0.0,
            }
        }
    }

    out_json.parent.mkdir(parents=True, exist_ok=True)
    out_json.write_text(json.dumps(summary, indent=2))
    print(f"✅ wrote {out_json}")

    if out_tsv:
        out_tsv.parent.mkdir(parents=True, exist_ok=True)
        out_df.to_csv(out_tsv, sep="\t", index=False)
        print(f"✅ wrote {out_tsv}")

def walk_tree_and_process(root: Path, contexts_list: List[str], dbnsfp: Path, out_name: str, write_tsv: bool):
    for sig_dir in sorted(root.glob("*")):
        if not sig_dir.is_dir():
            continue
        for n_dir in sorted(sig_dir.glob("n_*")):
            if not n_dir.is_dir():
                continue
            for parq in sorted(n_dir.glob("rep_*.sim.parquet")):
                m = re.search(r"rep_(\d+)", parq.stem)
                rep = f"rep_{m.group(1)}" if m else parq.stem
                out_json = n_dir / f"{out_name}_{rep}.json"
                out_tsv = n_dir / f"{out_name}_{rep}.tsv" if write_tsv else None
                process_replicate(parq, contexts_list, dbnsfp, out_json, out_tsv)

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--root", required=True, help="Root like results_DBS (contains DBS*/n_*/rep_*.sim.parquet)")
    ap.add_argument("--contexts_list", required=True, help="File listing contexts in context_id order (e.g. context_signature_list_DBS.txt)")
    ap.add_argument("--dbnsfp", required=True, help="bgzip+tabix dbNSFP file (e.g. ref/dbNSFP5.2a_grch38.gz)")
    ap.add_argument("--out-prefix", default="pathogenicity_summary", help="Prefix for output files (per replicate)")
    ap.add_argument("--write-tsv", action="store_true", help="Also write per-variant TSV per replicate (debug/QC)")
    args = ap.parse_args()

    root = Path(args.root)
    if not root.exists():
        raise SystemExit(f"❌ root not found: {root}")

    ctx_list = load_contexts_list(Path(args.contexts_list))
    if not ctx_list:
        raise SystemExit("❌ contexts_list is empty")

    dbnsfp = Path(args.dbnsfp)
    if not dbnsfp.exists():
        raise SystemExit(f"❌ dbNSFP not found: {dbnsfp}")

    walk_tree_and_process(root, ctx_list, dbnsfp, args.out_prefix, args.write_tsv)

if __name__ == "__main__":
    main()
