#!/usr/bin/env python3
from __future__ import annotations
import argparse, json, re, subprocess
from pathlib import Path

import pandas as pd
import pysam

# ================= thresholds (NUMERIC ONLY) =================
SIFT4G_DAMAGING_LT = 0.05                 # <0.05 damaging, else benign

PP2_BENIGN_LT      = 0.15                 # PolyPhen-2 HDIV
PP2_POSS_UPPER     = 0.85                 # [0.15, 0.85] possibly damaging
# >=0.85 damaging

CADD_BENIGN_LT     = 2.0                  # PHRED
CADD_UNCERT_LT     = 10.0                 # [2,10)
CADD_LIKELY_LT     = 20.0                 # [10,20)
# >=20 damaging

# ================= DBS78 list =================
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

def split_ctx(s: str) -> tuple[str,str]:
    a, b = s.split(">")
    return a, b

DNA_COMP = str.maketrans("ACGT","TGCA")
def revcomp(s: str) -> str:
    return s.translate(DNA_COMP)[::-1]

# ================= translation & consequence =================
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

# ================= dbNSFP (numeric only) =================
class DBNSFP:
    IDX_SIFT4G_SCORE   = 50 - 1
    IDX_PP2_HDIV_SCORE = 53 - 1
    IDX_CADD_PHRED     = 144 - 1

    def __init__(self, path: Path):
        self.path = str(path)

    def _split(self, s: str) -> list[str] | None:
        return s.split(';') if s and s != '.' else None

    def _to_float(self, s: str|None) -> float|None:
        try:
            return float(s) if s is not None else None
        except ValueError:
            return None

    def _median_from_list(self, vals: list[str]|None) -> float|None:
        if not vals:
            return None
        nums = []
        for v in vals:
            try: nums.append(float(v))
            except: pass
        if not nums:
            return None
        nums.sort()
        m = len(nums)//2
        return nums[m] if len(nums)%2 else (nums[m-1]+nums[m])/2

    def query(self, chrom: str, pos: int, ref: str, alt: str) -> dict:
        for chrom_try in (chrom.lstrip("chr"), f"chr{chrom.lstrip('chr')}"):
            try:
                res = subprocess.run(
                    ["tabix", self.path, f"{chrom_try}:{pos}-{pos}"],
                    capture_output=True, check=True
                )
            except subprocess.CalledProcessError:
                continue
            lines = [ln for ln in res.stdout.decode("utf-8","replace").splitlines() if ln]
            for ln in lines:
                f = ln.split("\t")
                if len(f) < 4: 
                    continue
                if f[2] != ref or f[3] != alt:
                    continue
                sift_scores = self._split(f[self.IDX_SIFT4G_SCORE])
                pp2_scores  = self._split(f[self.IDX_PP2_HDIV_SCORE])
                cadd_val    = self._to_float(f[self.IDX_CADD_PHRED])
                return {
                    "sift4g_score": self._median_from_list(sift_scores),
                    "pp2_hdiv_score": self._median_from_list(pp2_scores),
                    "cadd_phred": cadd_val
                }
        return {"sift4g_score": None, "pp2_hdiv_score": None, "cadd_phred": None}

# ----- numeric category helpers -----
def sift_cat(score: float|None) -> str:
    if score is None: return "missing"
    return "damaging" if score < SIFT4G_DAMAGING_LT else "benign"

def pp2_cat(score: float|None) -> str:
    if score is None: return "missing"
    if score < PP2_BENIGN_LT:  return "benign"
    if score < PP2_POSS_UPPER: return "possibly_damaging"
    return "damaging"

def cadd_cat(score: float|None) -> str:
    if score is None: return "missing"
    if score < CADD_BENIGN_LT:  return "benign"
    if score < CADD_UNCERT_LT:  return "uncertain"
    if score < CADD_LIKELY_LT:  return "likely_damaging"
    return "damaging"

# ================= AA consequence & SNV derivation =================
def aa_consequences_for_row(row, fasta: pysam.FastaFile) -> tuple[str,str|None,str]:
    chrom = str(row["chr"]); pos = int(row["pos"]); strand = int(row["strand"])
    ref_codon = str(row["ref_codon"]).upper(); ref_aa = str(row["ref_aa"]).upper()
    ci = int(row["codon_index"]); ctx = ID2CTX[int(row["context_id"])]
    _, alt2_ctx = split_ctx(ctx)

    if ci in (0,1):  # within a codon
        alt = list(ref_codon)
        alt[ci]   = alt2_ctx[0]
        alt[ci+1] = alt2_ctx[1]
        lbl = classify(ref_aa, aa_of("".join(alt)))
        return lbl, None, "left"

    # boundary: need next codon on transcript strand
    if strand == 1:
        next_start = pos + 1
        next_trip  = fasta.fetch(chrom, next_start-1, next_start-1+3).upper()
    else:
        next_start = pos - 1
        seq = fasta.fetch(chrom, next_start-3, next_start).upper()
        next_trip = revcomp(seq)

    alt_left  = list(ref_codon); alt_left[2] = alt2_ctx[0]; alt_left = "".join(alt_left)
    alt_right = list(next_trip);  alt_right[0] = alt2_ctx[1]; alt_right = "".join(alt_right)

    lbl_left  = classify(ref_aa, aa_of(alt_left))
    lbl_right = classify(aa_of(next_trip) or "X", aa_of(alt_right))
    primary   = worse(lbl_left, lbl_right)
    secondary = lbl_right if primary == lbl_left else lbl_left
    primary_side = "left" if primary == lbl_left else "right"
    return primary, secondary, primary_side

def forward_refalt_for_row(row, fasta: pysam.FastaFile):
    chrom = str(row["chr"]); pos = int(row["pos"])
    gref2 = fasta.fetch(chrom, pos-1, pos+1).upper()
    ctx = ID2CTX[int(row["context_id"])]
    ref2_ctx, alt2_ctx = split_ctx(ctx)
    if ref2_ctx != gref2:
        if revcomp(ref2_ctx) == gref2:
            ref2_ctx = gref2
            alt2_ctx = revcomp(alt2_ctx)
        else:
            alt2_ctx = gref2
    snv1 = (chrom, pos,     gref2[0], alt2_ctx[0])
    snv2 = (chrom, pos + 1, gref2[1], alt2_ctx[1])
    return snv1, snv2

def aggregate_same_codon(s1: dict, s2: dict) -> dict:
    def max_or_none(a,b):
        vals=[v for v in (a,b) if v is not None]
        return max(vals) if vals else None
    def min_or_none(a,b):
        vals=[v for v in (a,b) if v is not None]
        return min(vals) if vals else None
    return {
        "sift4g_score":  min_or_none(s1["sift4g_score"],  s2["sift4g_score"]),   # lower = worse
        "pp2_hdiv_score":max_or_none(s1["pp2_hdiv_score"],s2["pp2_hdiv_score"]), # higher = worse
        "cadd_phred":    max_or_none(s1["cadd_phred"],    s2["cadd_phred"])      # higher = worse
    }

# ================= summarization =================
def summarize_replicate(pq: Path, fasta: pysam.FastaFile, db: "DBNSFP", write_tsv: bool):
    rep_id = pq.stem
    cols = ["chr","pos","strand","ref_codon","ref_aa","codon_index","context_id"]
    df = pd.read_parquet(pq, columns=cols)
    if df.empty:
        return

    # category counters
    sift_counts  = {"damaging":0, "benign":0, "missing":0}
    pp2_counts   = {"benign":0, "possibly_damaging":0, "damaging":0, "missing":0}
    cadd_counts  = {"benign":0, "uncertain":0, "likely_damaging":0, "damaging":0, "missing":0}
    n = 0

    out_rows = []

    for _, row in df.iterrows():
        primary_lbl, secondary_lbl, primary_side = aa_consequences_for_row(row, fasta)
        snv1, snv2 = forward_refalt_for_row(row, fasta)
        s1 = db.query(*snv1)
        s2 = db.query(*snv2)

        if int(row["codon_index"]) in (0,1):
            sc = aggregate_same_codon(s1, s2)
        else:
            sc = s1 if primary_side == "left" else s2

        n += 1
        sift_counts[ sift_cat(sc["sift4g_score"]) ] += 1
        pp2_counts[  pp2_cat(sc["pp2_hdiv_score"]) ] += 1
        cadd_counts[ cadd_cat(sc["cadd_phred"]) ]    += 1

        if write_tsv:
            out_rows.append({
                "chr": row["chr"], "pos": int(row["pos"]), "codon_index": int(row["codon_index"]),
                "primary_label": primary_lbl, "secondary_label": secondary_lbl, "primary_side": primary_side,
                "sift4g_score": sc["sift4g_score"],
                "pp2_hdiv_score": sc["pp2_hdiv_score"],
                "cadd_phred": sc["cadd_phred"]
            })

    summary = {
        "replicate": rep_id,
        "n_variants": n,
        "SIFT4G": {
            "counts": sift_counts,
            "fractions": {k: (v/max(n,1)) for k,v in sift_counts.items()},
            "thresholds": {"damaging": f"< {SIFT4G_DAMAGING_LT}", "benign": f">= {SIFT4G_DAMAGING_LT}"}
        },
        "PolyPhen2_HDIV": {
            "counts": pp2_counts,
            "fractions": {k: (v/max(n,1)) for k,v in pp2_counts.items()},
            "thresholds": {
                "benign": f"< {PP2_BENIGN_LT}",
                "possibly_damaging": f"[{PP2_BENIGN_LT}, {PP2_POSS_UPPER}]",
                "damaging": f">= {PP2_POSS_UPPER}"
            }
        },
        "CADD_PHRED": {
            "counts": cadd_counts,
            "fractions": {k: (v/max(n,1)) for k,v in cadd_counts.items()},
            "thresholds": {
                "benign": f"< {CADD_BENIGN_LT}",
                "uncertain": f"[{CADD_BENIGN_LT}, {CADD_UNCERT_LT})",
                "likely_damaging": f"[{CADD_UNCERT_LT}, {CADD_LIKELY_LT})",
                "damaging": f">= {CADD_LIKELY_LT}"
            }
        }
    }

    out_json = pq.parent / f"pathogenicity_summary_{rep_id}.json"
    out_json.write_text(json.dumps(summary, indent=2))

    if write_tsv and out_rows:
        out_tsv = pq.parent / f"pathogenicity_primary_scores_{rep_id}.tsv"
        pd.DataFrame(out_rows).to_csv(out_tsv, sep="\t", index=False)

    print(f"✅ {pq.parent.name}/{rep_id}: wrote {out_json.name}")

# ================= main =================
def main():
    ap = argparse.ArgumentParser(description="Compile NUMERIC dbNSFP scores for DBS sims (AA-centric).")
    ap.add_argument("--root", required=True, help="Root with DBS*/n_*/rep_*.sim.with_aapos.parquet")
    ap.add_argument("--fasta", required=True, help="Reference FASTA (indexed .fai)")
    ap.add_argument("--dbnsfp", required=True, help="bgzip+tabix dbNSFP file")
    ap.add_argument("--write-tsv", action="store_true", help="Also write per-variant TSV with primary numeric scores")
    args = ap.parse_args()

    root  = Path(args.root)
    fasta = pysam.FastaFile(args.fasta)
    db    = DBNSFP(Path(args.dbnsfp))

    for sig_dir in sorted(root.glob("DBS*")):
        for n_dir in sorted(sig_dir.glob("n_*")):
            for pq in sorted(n_dir.glob("rep_*.sim.with_aapos.parquet")):
                summarize_replicate(pq, fasta, db, args.write_tsv)

if __name__ == "__main__":
    main()
