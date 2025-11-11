#!/usr/bin/env python3
"""
Annotate simulated DBS mutations at the amino acid level without VEP.

For each replicate (rep_*.sim.parquet):
- Determine if each DBS changes one or two codons.
- Translate altered codons to AA(s).
- Assign the most deleterious consequence:
    stop_gained > missense > synonymous
- If two codons are affected, record both (primary=worst, secondary=other).
- Write one JSON summary per replicate: consequence_summary_rep_X.json
- Optional per-variant TSV with primary/secondary calls.

Usage:
  python add_csq_v2.py \
      --root results_DBS \
      --fasta ref/GRCh38.primary_assembly.genome.fa \
      [--write-tsv]
"""

from __future__ import annotations
import argparse, json, re
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

# ---------------- DNA helpers ----------------
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

# ---------------- per-variant call ----------------
def call_consequence_row(row, fasta: pysam.FastaFile) -> tuple[str,str|None]:
    chrom = str(row["chr"])
    pos   = int(row["pos"])
    strand = int(row["strand"])
    ref_codon = str(row["ref_codon"]).upper()
    ref_aa     = str(row["ref_aa"]).upper()
    ci         = int(row["codon_index"])
    ctx        = ID2CTX[int(row["context_id"])]
    ref2, alt2 = split_ctx(ctx)

    if strand == -1:
        ref2 = revcomp(ref2)
        alt2 = revcomp(alt2)

    # Single-codon cases
    if ci in (0,1):
        alt_codon = list(ref_codon)
        alt_codon[ci]   = alt2[0]
        alt_codon[ci+1] = alt2[1]
        alt_codon = "".join(alt_codon)
        alt_aa = aa_of(alt_codon)
        return classify(ref_aa, alt_aa or ref_aa), None

    # Boundary case (ci == 2)
    if strand == 1:
        next_start = pos + 1
        next_trip = fasta.fetch(chrom, next_start-1, next_start-1+3).upper()
    else:
        next_start = pos - 1
        seq = fasta.fetch(chrom, next_start-3, next_start).upper()
        next_trip = revcomp(seq)

    if len(next_trip) != 3 or "N" in next_trip:
        alt_codon = list(ref_codon)
        alt_codon[2] = alt2[0]
        alt_codon = "".join(alt_codon)
        alt_aa = aa_of(alt_codon)
        return classify(ref_aa, alt_aa or ref_aa), None

    alt1 = list(ref_codon); alt1[2] = alt2[0]; alt1 = "".join(alt1)
    alt2cod = list(next_trip); alt2cod[0] = alt2[1]; alt2cod = "".join(alt2cod)

    aa1 = aa_of(alt1) or ref_aa
    aa2 = aa_of(alt2cod) or aa_of(next_trip) or "X"
    label1 = classify(ref_aa, aa1)
    label2 = classify(aa_of(next_trip) or "X", aa2)
    primary, secondary = worst(label1, label2)
    return primary, secondary

# ---------------- summarize per replicate ----------------
def summarize_replicate(pq: Path, fasta: pysam.FastaFile, write_tsv: bool):
    rep_id = pq.stem
    cols = ["chr","pos","strand","ref_codon","ref_aa","codon_index","context_id"]
    df = pd.read_parquet(pq, columns=cols)
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
                "codon_index": int(row["codon_index"]),
                "ref_codon": row["ref_codon"], "ref_aa": row["ref_aa"],
                "context_id": int(row["context_id"]),
                "primary": primary, "secondary": secondary
            })

    n_total = len(df)
    summary = {
        "replicate": rep_id,
        "n_mutations": n_total,
        "counts": dict(c_primary),
        "fractions": {k: round(v/n_total,6) for k,v in c_primary.items()},
        "secondary_counts": dict(c_secondary),
    }

    out_json = pq.parent / f"consequence_summary_{rep_id}.json"
    out_json.write_text(json.dumps(summary, indent=2))

    if write_tsv and per_var:
        out_tsv = pq.parent / f"consequence_per_variant_{rep_id}.tsv"
        pd.DataFrame(per_var).to_csv(out_tsv, sep="\t", index=False)

    print(f"✅ {rep_id}: {n_total} variants → {out_json.name}")
    return summary

# ---------------- main ----------------
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
            for pq in sorted(n_dir.glob("rep_*.sim.parquet")):
                summarize_replicate(pq, fasta, args.write_tsv)

if __name__ == "__main__":
    main()
