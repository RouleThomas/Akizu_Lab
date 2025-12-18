#!/usr/bin/env python3
"""
Build a per-nucleotide CDS (coding) table for one chromosome using:
  - GTF (CDS features with transcript structure)
  - Reference genome FASTA

Outputs Parquet with 1-based genomic positions (IGV-friendly).

Columns:
  chr, pos(1-based), strand(+1/-1), ref_base(0:A 1:C 2:G 3:T 4:N),
  gene_id, transcript_id,
  cds_offset(0-based in spliced CDS),
  codon_index(0-based codon number),
  codon_pos(0/1/2 within codon),
  ref_codon(3-mer in coding orientation),
  ref_aa(single letter; '*' for stop; 'X' for unknown)

This FIXES:
  - splicing (exon/CDS structure)
  - minus-strand orientation
  - "frame shifts" caused by resetting i%3 per exon
"""

from __future__ import annotations
import argparse, gzip, re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple, Iterable, Optional

import pyarrow as pa
import pyarrow.parquet as pq

try:
    import pysam
except ImportError as e:
    raise SystemExit("ERROR: pysam is required. Install it in your env (pip/conda).") from e

from Bio.Data import CodonTable

# --- constants ---
CODON_TABLE = CodonTable.unambiguous_dna_by_name["Standard"].forward_table
STOP_CODONS = {"TAA", "TAG", "TGA"}
BASE2INT = {"A": 0, "C": 1, "G": 2, "T": 3, "N": 4}

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_attrs(attr: str) -> Dict[str, str]:
    return {k: v for k, v in ATTR_RE.findall(attr)}

def revcomp(seq: str) -> str:
    comp = str.maketrans("ACGTNacgtn", "TGCANtgcan")
    return seq.translate(comp)[::-1]

@dataclass
class CDSExon:
    start: int  # 1-based inclusive
    end: int    # 1-based inclusive
    phase: int  # 0/1/2 or -1

def open_text_maybe_gz(path: str):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path, "rt")

def aa_for_codon(c: str) -> str:
    c = c.upper()
    if "N" in c or len(c) != 3:
        return "X"
    if c in STOP_CODONS:
        return "*"
    return CODON_TABLE.get(c, "X")

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--gtf", required=True, help="GTF or GTF.GZ with CDS features (e.g., gencode)")
    ap.add_argument("--genome", required=True, help="Reference genome FASTA (must be faidx-indexed)")
    ap.add_argument("--chrom", required=True, help="Chromosome name (e.g., chr1)")
    ap.add_argument("--output", required=True, help="Output parquet path")
    ap.add_argument("--flush_every", type=int, default=2_000_000, help="Rows per parquet row-group flush")
    args = ap.parse_args()

    out_path = Path(args.output)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    # Load genome
    fa = pysam.FastaFile(args.genome)

    # Collect CDS intervals per transcript on this chrom
    # transcript_id -> (gene_id, strand, [CDSExon...])
    cds_by_tx: Dict[str, Tuple[str, int, List[CDSExon]]] = {}

    tmp_exons: Dict[str, List[CDSExon]] = defaultdict(list)
    tx_gene: Dict[str, str] = {}
    tx_strand: Dict[str, int] = {}

    with open_text_maybe_gz(args.gtf) as fh:
        for ln in fh:
            if not ln or ln[0] == "#":
                continue
            parts = ln.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, source, feature, start_s, end_s, score, strand_s, frame_s, attrs = parts
            if chrom != args.chrom:
                continue
            if feature != "CDS":
                continue

            a = parse_attrs(attrs)
            tid = a.get("transcript_id")
            gid = a.get("gene_id")
            if not tid or not gid:
                continue

            start = int(start_s)
            end = int(end_s)
            strand = 1 if strand_s == "+" else -1
            phase = int(frame_s) if frame_s.isdigit() else -1

            tmp_exons[tid].append(CDSExon(start=start, end=end, phase=phase))
            tx_gene[tid] = gid
            tx_strand[tid] = strand

    # Freeze into cds_by_tx
    for tid, exons in tmp_exons.items():
        cds_by_tx[tid] = (tx_gene[tid], tx_strand[tid], exons)

    # Parquet writer (streaming)
    schema = pa.schema([
        ("chr", pa.string()),
        ("pos", pa.int32()),              # 1-based genomic
        ("strand", pa.int8()),            # +1/-1
        ("ref_base", pa.int8()),          # A/C/G/T/N -> 0/1/2/3/4
        ("gene_id", pa.string()),
        ("transcript_id", pa.string()),
        ("cds_offset", pa.int32()),       # 0-based in spliced CDS
        ("codon_index", pa.int32()),      # 0-based codon number
        ("codon_pos", pa.int8()),         # 0/1/2 within codon
        ("ref_codon", pa.string()),       # 3-mer in coding orientation
        ("ref_aa", pa.string()),          # single-letter
    ])

    writer = pq.ParquetWriter(out_path, schema=schema, compression="zstd")

    # Buffers
    buf = {name: [] for name in schema.names}
    n_buf = 0

    def flush():
        nonlocal n_buf
        if n_buf == 0:
            return
        table = pa.Table.from_pydict(buf, schema=schema)
        writer.write_table(table)
        for k in buf:
            buf[k].clear()
        n_buf = 0

    # Build per transcript: spliced CDS sequence + mapping
    for tid, (gid, strand, exons) in cds_by_tx.items():
        # sort exons in transcript order
        exons_sorted = sorted(exons, key=lambda e: e.start, reverse=(strand == -1))

        cds_seq_parts: List[str] = []
        mapping: List[Tuple[int, str]] = []  # list of (genomic_pos_1based, base_in_coding_orientation)
        # Build spliced CDS in coding orientation (5'->3' of transcript CDS)
        for ex in exons_sorted:
            # pysam fetch is 0-based half-open, so convert: [start-1, end)
            seg = fa.fetch(args.chrom, ex.start - 1, ex.end).upper()
            if strand == -1:
                seg = revcomp(seg)
                # transcript direction goes from high->low genomic for '-' :
                # mapping positions should follow coding orientation
                pos_iter = range(ex.end, ex.start - 1, -1)  # end..start
            else:
                pos_iter = range(ex.start, ex.end + 1)      # start..end

            # Add bases in transcript/coding order
            cds_seq_parts.append(seg)
            for p, b in zip(pos_iter, seg):
                mapping.append((p, b))

        cds_seq = "".join(cds_seq_parts)
        if not cds_seq:
            continue

        # (Optional sanity) verify phase consistency if available (won't stop, just helps debugging)
        # Phase meaning is tricky; we won't enforce it hard. The crucial part is: we DO NOT reset frame per exon.

        # Emit rows
        for cds_offset, (pos1, base) in enumerate(mapping):
            codon_pos = cds_offset % 3
            codon_start = cds_offset - codon_pos
            codon = cds_seq[codon_start:codon_start + 3]
            if len(codon) < 3:
                codon = codon + ("N" * (3 - len(codon)))

            buf["chr"].append(args.chrom)
            buf["pos"].append(int(pos1))  # 1-based
            buf["strand"].append(int(strand))
            buf["ref_base"].append(int(BASE2INT.get(base, 4)))
            buf["gene_id"].append(gid)
            buf["transcript_id"].append(tid)
            buf["cds_offset"].append(int(cds_offset))
            buf["codon_index"].append(int(codon_start // 3))
            buf["codon_pos"].append(int(codon_pos))
            buf["ref_codon"].append(codon)
            buf["ref_aa"].append(aa_for_codon(codon))

            n_buf += 1
            if n_buf >= args.flush_every:
                flush()

    flush()
    writer.close()
    print(f"✅ wrote {out_path}")

if __name__ == "__main__":
    main()
