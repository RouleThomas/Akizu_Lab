#!/usr/bin/env python3
"""
Annotate simulated DBS mutations at the amino acid level without VEP.

Assumptions (matches your new DBS pipeline):
- Parquet has 1-based genomic pos.
- DBS is anchored at genomic (pos, pos+1) on the forward genome.
- context_name (e.g. "AC>CA") is forward-genome representation.
- strand is transcript strand (+1 / -1).
- codon_pos is 0/1/2 position of THIS row's base within ref_codon in transcript orientation.
  (If codon_pos missing, we fall back to codon_index modulo 3.)
"""

from __future__ import annotations
import argparse, json
from pathlib import Path
from collections import Counter

import pandas as pd
import pysam

# ---------------- DBS contexts ----------------
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
    ref2, alt2 = s.split(">")
    return ref2, alt2

DNA_COMP = str.maketrans("ACGT", "TGCA")
def revcomp(s: str) -> str:
    return s.translate(DNA_COMP)[::-1]

# ---------------- translation table ----------------
CODON2AA = {
    "TTT":"F","TTC":"F","TTA":"L","TTG":"L",
    "CTT":"L","CTC":"L","CTA":"L","CTG":"L",
    "ATT":"I","ATC":"I","ATA":"I","ATG":"M",
    "GTT":"V","GTC":"V","GTA":"V","GTG":"V",
    "TCT":"S","TCC":"S","TCA":"S","TCG":"S",
    "CCT":"P","CCC":"P","CCA":"P","CCG":"P",
    "ACT":"T","ACC":"T","ACA":"T","ACG":"T",
    "GCT":"A","GCC":"A","GCA":"A","GCG":"A",
    "TAT":"Y","TAC":"Y","TAA":"*","TAG":"*",
    "CAT":"H","CAC":"H","CAA":"Q","CAG":"Q",
    "AAT":"N","AAC":"N","AAA":"K","AAG":"K",
    "GAT":"D","GAC":"D","GAA":"E","GAG":"E",
    "TGT":"C","TGC":"C","TGA":"*","TGG":"W",
    "CGT":"R","CGC":"R","CGA":"R","CGG":"R",
    "AGT":"S","AGC":"S","AGA":"R","AGG":"R",
    "GGT":"G","GGC":"G","GGA":"G","GGG":"G"
}

def aa_of(codon: str) -> str|None:
    codon = codon.upper().replace("U","T")
    return CODON2AA.get(codon)

def classify(ref_aa: str, alt_aa: str) -> str:
    if alt_aa == "*":
        return "stop_gained"
    if alt_aa != ref_aa:
        return "missense"
    return "synonymous"

def worst(primary: str, secondary: str|None) -> tuple[str,str|None]:
    order = {"stop_gained":0,"missense":1,"synonymous":2}
    if secondary is None:
        return primary, None
    return (primary, secondary) if order[primary] <= order[secondary] else (secondary, primary)

# ---- fetch codon on transcript orientation using genome ----
def fetch_triplet_plus(fasta, chrom: str, start_1based: int) -> str:
    # + strand codon at start_1based
    return fasta.fetch(chrom, start_1based-1, start_1based-1+3).upper()

def fetch_triplet_minus(fasta, chrom: str, start_1based: int) -> str:
    # transcript codon starting at start_1based on transcript corresponds to genomic [start-2 .. start] (inclusive)
    # easiest: fetch genomic (start_1based-2 .. start_1based) then revcomp
    seq = fasta.fetch(chrom, start_1based-3, start_1based).upper()
    return revcomp(seq)

# ---------------- per-variant call ----------------
def call_consequence_row(row, fasta: pysam.FastaFile) -> tuple[str,str|None]:
    chrom = str(row["chr"])
    pos   = int(row["pos"])      # 1-based
    strand = int(row["strand"])
    ref_codon = str(row["ref_codon"]).upper()
    ref_aa     = str(row["ref_aa"]).upper()

    # Prefer codon_pos (your new parquet). Fallback: codon_index % 3.
    if "codon_pos" in row and pd.notna(row["codon_pos"]):
        cp = int(row["codon_pos"])
    else:
        cp = int(row["codon_index"]) % 3

    # Context string: prefer context_name if present
    if "context_name" in row and pd.notna(row["context_name"]):
        ctx = str(row["context_name"])
    else:
        ctx = ID2CTX[int(row["context_id"])]

    ref2, alt2 = split_ctx(ctx)

    # Convert forward-genome representation into transcript-orientation representation
    # For strand -1, the dinucleotide in transcript order is revcomp(forward_genome_dinuc)
    if strand == -1:
        ref2 = revcomp(ref2)
        alt2 = revcomp(alt2)

    # Now alt2 is in transcript order relative to the affected transcript bases.
    # BUT: our row is anchored at genomic pos (the left base of forward genome).
    # Mapping of which alt base hits THIS row depends on strand:
    #
    # strand +1: this row corresponds to the FIRST base of the dinuc in transcript order -> alt2[0]
    # strand -1: this row corresponds to the SECOND base in transcript order -> alt2[1]
    if strand == 1:
        alt_for_this = alt2[0]
        alt_for_adj  = alt2[1]  # adjacent transcript base is cp+1 (or next codon)
    else:
        alt_for_this = alt2[1]
        alt_for_adj  = alt2[0]  # adjacent transcript base is cp-1 (or previous codon)

    # Basic guard
    if len(ref_codon) != 3 or "N" in ref_codon:
        return "synonymous", None

    # ---------------- Single-codon cases ----------------
    # strand +1 affects (cp, cp+1) within same codon if cp in {0,1}
    if strand == 1 and cp in (0, 1):
        alt_codon = list(ref_codon)
        alt_codon[cp]   = alt_for_this
        alt_codon[cp+1] = alt_for_adj
        alt_codon = "".join(alt_codon)
        alt_aa = aa_of(alt_codon) or ref_aa
        return classify(ref_aa, alt_aa), None

    # strand -1 affects (cp-1, cp) within same codon if cp in {1,2}
    if strand == -1 and cp in (1, 2):
        alt_codon = list(ref_codon)
        alt_codon[cp]   = alt_for_this
        alt_codon[cp-1] = alt_for_adj
        alt_codon = "".join(alt_codon)
        alt_aa = aa_of(alt_codon) or ref_aa
        return classify(ref_aa, alt_aa), None

    # ---------------- Boundary cases ----------------
    # Case A: strand +1 and cp == 2 => affects last base of this codon + first base of NEXT codon
    if strand == 1 and cp == 2:
        # mutate this codon at pos 2 with alt_for_this
        alt1 = list(ref_codon); alt1[2] = alt_for_this; alt1 = "".join(alt1)

        # next codon on transcript starts at genomic pos+1
        next_trip = fetch_triplet_plus(fasta, chrom, pos + 1)
        if len(next_trip) != 3 or "N" in next_trip:
            aa1 = aa_of(alt1) or ref_aa
            return classify(ref_aa, aa1), None

        alt2cod = list(next_trip); alt2cod[0] = alt_for_adj; alt2cod = "".join(alt2cod)

        aa1_ref = ref_aa
        aa1_alt = aa_of(alt1) or aa1_ref

        aa2_ref = aa_of(next_trip) or "X"
        aa2_alt = aa_of(alt2cod) or aa2_ref

        label1 = classify(aa1_ref, aa1_alt)
        label2 = classify(aa2_ref, aa2_alt)
        return worst(label1, label2)

    # Case B: strand -1 and cp == 0 => affects first base of this codon + last base of PREVIOUS codon
    if strand == -1 and cp == 0:
        # mutate this codon at pos 0 with alt_for_this
        alt1 = list(ref_codon); alt1[0] = alt_for_this; alt1 = "".join(alt1)

        # previous codon on transcript corresponds to genomic pos+1 as the codon "start" in transcript direction
        # Its transcript triplet can be obtained by fetching genomic [pos+1 .. pos+3] then revcomp.
        prev_trip = fasta.fetch(chrom, (pos+1)-1, (pos+1)-1+3).upper()
        prev_trip = revcomp(prev_trip)

        if len(prev_trip) != 3 or "N" in prev_trip:
            aa1 = aa_of(alt1) or ref_aa
            return classify(ref_aa, aa1), None

        alt2cod = list(prev_trip); alt2cod[2] = alt_for_adj; alt2cod = "".join(alt2cod)

        aa1_ref = ref_aa
        aa1_alt = aa_of(alt1) or aa1_ref

        aa2_ref = aa_of(prev_trip) or "X"
        aa2_alt = aa_of(alt2cod) or aa2_ref

        label1 = classify(aa1_ref, aa1_alt)
        label2 = classify(aa2_ref, aa2_alt)
        return worst(label1, label2)

    # If we fall through, be conservative
    return "synonymous", None


# ---------------- summarize per replicate ----------------
def summarize_replicate(pq_path: Path, fasta: pysam.FastaFile, write_tsv: bool):
    rep_id = pq_path.stem

    cols = ["chr","pos","strand","ref_codon","ref_aa","context_id"]
    # new parquet usually has codon_pos + context_name; keep optional
    opt_cols = ["codon_pos","codon_index","context_name"]
    df = pd.read_parquet(pq_path, columns=cols + opt_cols)

    # coding-only
    coding = df["ref_codon"].fillna("").astype(str).str.len() == 3
    df = df.loc[coding].copy()
    if df.empty:
        return None

    per_var = []
    c_primary, c_secondary = Counter(), Counter()

    for _, row in df.iterrows():
        primary, secondary = call_consequence_row(row, fasta)
        c_primary[primary] += 1
        if secondary:
            c_secondary[secondary] += 1

        if write_tsv:
            per_var.append({
                "chr": row["chr"], "pos": int(row["pos"]), "strand": int(row["strand"]),
                "codon_pos": (int(row["codon_pos"]) if pd.notna(row.get("codon_pos", None)) else None),
                "codon_index": (int(row["codon_index"]) if pd.notna(row.get("codon_index", None)) else None),
                "ref_codon": row["ref_codon"], "ref_aa": row["ref_aa"],
                "context_id": int(row["context_id"]),
                "context_name": (str(row["context_name"]) if pd.notna(row.get("context_name", None)) else None),
                "primary": primary, "secondary": secondary
            })

    n_total = len(df)
    summary = {
        "replicate": rep_id,
        "n_mutations": int(n_total),
        "counts_primary": dict(c_primary),
        "fractions_primary": {k: round(v/n_total, 6) for k, v in c_primary.items()},
        "counts_secondary": dict(c_secondary),
    }

    out_json = pq_path.parent / f"consequence_summary_{rep_id}.json"
    out_json.write_text(json.dumps(summary, indent=2))

    if write_tsv and per_var:
        out_tsv = pq_path.parent / f"consequence_per_variant_{rep_id}.tsv"
        pd.DataFrame(per_var).to_csv(out_tsv, sep="\t", index=False)

    print(f"✅ {rep_id}: {n_total} coding DBS → {out_json.name}")
    return summary


def main():
    ap = argparse.ArgumentParser(description="AA-level DBS consequence caller per replicate")
    ap.add_argument("--root", required=True, help="Root folder (DBS*/n_*/rep_*.sim.parquet)")
    ap.add_argument("--fasta", required=True, help="Reference genomic FASTA (indexed)")
    ap.add_argument("--write-tsv", action="store_true", help="Also write per-variant TSVs")
    args = ap.parse_args()

    fasta = pysam.FastaFile(args.fasta)
    root = Path(args.root)

    for sig_dir in sorted(root.glob("*")):
        for n_dir in sorted(sig_dir.glob("n_*")):
            for pq_path in sorted(n_dir.glob("rep_*.sim.parquet")):
                summarize_replicate(pq_path, fasta, args.write_tsv)

if __name__ == "__main__":
    main()
