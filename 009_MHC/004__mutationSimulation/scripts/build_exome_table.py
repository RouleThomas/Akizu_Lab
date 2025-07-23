#!/usr/bin/env python3
"""
Create a per-nucleotide Parquet file for one chromosome from a CDS FASTA file.
Each FASTA header must be in this format:
    >ENSG00000186092.7::chr1:65564-65573(+)

This script outputs one compressed .parquet file per chromosome,
always using the + strand genomic coordinates (no reverse complement).
"""

from __future__ import annotations
import argparse, re
from pathlib import Path
from typing import Dict, List

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO
from Bio.Data import CodonTable

# ------------------------------------------------------------------
# Constants
# ------------------------------------------------------------------
CODON_TABLE: Dict[str, str] = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
STOP_CODONS = {"TAA", "TAG", "TGA"}
BASE2INT = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}
HEADER_RE = re.compile(r"^(ENSG[^:]+)::(chr[\w]+):(\d+)-(\d+)\(([+-])\)")

# ------------------------------------------------------------------
def parse_args():
    p = argparse.ArgumentParser()
    p.add_argument("--cds-fa", required=True)
    p.add_argument("--chrom", required=True)
    p.add_argument("--output", required=True)
    return p.parse_args()

# ------------------------------------------------------------------
def emit_row(chr_: str, pos: int, strand: int, base: str, gene: str,
             codon_idx: int, codon: str) -> Dict:
    aa = CODON_TABLE.get(codon, "*" if codon in STOP_CODONS else "N")
    return {
        "chr": chr_,
        "pos": int(pos),
        "strand": int(strand),
        "ref_base": int(BASE2INT.get(base, 4)),  # default to 'N'
        "gene_id": gene,
        "transcript_id": None,
        "codon_index": int(codon_idx),
        "ref_codon": codon,
        "ref_aa": aa,
    }

# ------------------------------------------------------------------
def main():
    args = parse_args()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: List[Dict] = []
    n_skipped = 0

    for rec in SeqIO.parse(args.cds_fa, "fasta"):
        m = HEADER_RE.match(rec.description)
        if not m:
            n_skipped += 1
            continue

        gene, chr_, start_s, end_s, strand_char = m.groups()
        if chr_ != args.chrom:
            continue

        start, end = int(start_s), int(end_s)
        strand = 1 if strand_char == "+" else -1

        seq = str(rec.seq).upper()

        for i, base in enumerate(seq):
            genomic_pos = start + i  # always increasing on + strand
            codon_idx = i % 3
            codon_start = i - codon_idx
            codon_seq = seq[codon_start: codon_start + 3]

            if len(codon_seq) < 3:
                codon_seq += "N" * (3 - len(codon_seq))

            rows.append(emit_row(chr_, genomic_pos, strand, base, gene,
                                 codon_idx, codon_seq))

    if not rows:
        raise ValueError(f"No CDS entries matched {args.chrom}")

    tbl = pa.Table.from_pylist(rows)
    pq.write_table(tbl, out_path, compression="zstd")
    print(f"✅ Wrote {tbl.num_rows:,} rows to {out_path}")
    if n_skipped > 0:
        print(f"⚠️  Skipped {n_skipped:,} CDS entries due to malformed headers")

# ------------------------------------------------------------------
if __name__ == "__main__":
    main()
