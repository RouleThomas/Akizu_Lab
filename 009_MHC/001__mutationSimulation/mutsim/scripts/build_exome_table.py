# build_exome_table.py (v3 – write once, no Parquet append)
"""
Create a per‑nucleotide Parquet file for one chromosome from a CDS FASTA whose
headers look like:
    >ENSG00000186092.7::chr1:65564-65573(+)

This version **loads the entire chromosome worth of CDS bases into memory and
writes a single Parquet file at the end** (simpler than chunk‑append, avoids
append issues with older PyArrow).

Usage (inside SLURM array):
    python build_exome_table.py \
       --cds-fa  ~/mutsim/gtf/cds.fa \
       --chrom   chr1 \
       --output  ~/mutsim/parquet/chr1.parquet

Estimated RAM:  <2 GB for the largest chromosome (chr1).
"""
from __future__ import annotations

import argparse
import re
from pathlib import Path
from typing import Dict, List

import numpy as np
import pyarrow as pa
import pyarrow.parquet as pq
from Bio import SeqIO
from Bio.Data import CodonTable

# ------------------------------------------------------------------
# constants / look‑ups
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
        "pos": pos,
        "strand": strand,
        "ref_base": BASE2INT.get(base, 4),
        "gene_id": gene,
        "transcript_id": None,
        "codon_index": codon_idx,
        "ref_codon": codon,
        "ref_aa": aa,
    }

# ------------------------------------------------------------------

def main():
    args = parse_args()
    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    rows: List[Dict] = []

    for rec in SeqIO.parse(args.cds_fa, "fasta"):
        m = HEADER_RE.match(rec.description)
        if not m:
            continue
        gene, chr_, start_s, end_s, strand_char = m.groups()
        if chr_ != args.chrom:
            continue
        start = int(start_s)
        end = int(end_s)
        strand = 1 if strand_char == "+" else -1

        seq = str(rec.seq).upper()
        if strand == -1:
            seq = str(rec.seq.reverse_complement()).upper()

        for i, base in enumerate(seq):
            genomic_pos = start + i if strand == 1 else end - i
            codon_idx = i % 3
            codon_start = i - codon_idx
            codon_seq = seq[codon_start: codon_start + 3]
            if len(codon_seq) < 3:
                codon_seq += "N" * (3 - len(codon_seq))
            rows.append(emit_row(chr_, genomic_pos, strand, base, gene,
                                   codon_idx, codon_seq))

    tbl = pa.Table.from_pylist(rows)
    pq.write_table(tbl, out_path, compression="zstd")
    print(f"wrote {tbl.num_rows:,} rows to {out_path}")


if __name__ == "__main__":
    main()

