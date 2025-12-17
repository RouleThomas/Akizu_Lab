#!/usr/bin/env python3
import argparse, gzip, pickle, re
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import pandas as pd

TX_ID_RE = re.compile(r'transcript_id "([^"]+)"')
GENE_ID_RE = re.compile(r'gene_id "([^"]+)"')

def parse_attr(attr: str, regex: re.Pattern) -> Optional[str]:
    m = regex.search(attr)
    return m.group(1) if m else None

def strip_version(x: str) -> str:
    # ENST00000.13 -> ENST00000
    return re.sub(r"\.\d+$", "", x) if x else x

def norm_chr(c: str) -> str:
    return c[3:] if c.startswith("chr") else c

@dataclass
class CDSModel:
    chrom: str          # normalized (no "chr")
    strand: int         # +1/-1
    exons: List[Tuple[int, int]]      # (start,end) 1-based inclusive, transcription order
    cumlen: List[int]                 # cds bases before each exon
    gene_id: str
    tx_id: str

def build_models_from_gtf(gtf_path: Path):
    opener = gzip.open if str(gtf_path).endswith(".gz") else open

    # collect CDS exons per transcript
    cds_by_tx: Dict[str, Dict[Tuple[str,int], List[Tuple[int,int]]]] = {}
    gene_by_tx: Dict[str, str] = {}

    with opener(gtf_path, "rt") as f:
        for line in f:
            if not line or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, feature, start, end, _, strand, _, attr = parts
            if feature != "CDS":
                continue

            tx = parse_attr(attr, TX_ID_RE)
            gene = parse_attr(attr, GENE_ID_RE)
            if not tx or not gene:
                continue

            tx0 = strip_version(tx)
            gene0 = strip_version(gene)

            st, en = int(start), int(end)
            sgn = 1 if strand == "+" else -1
            chrom_n = norm_chr(chrom)

            cds_by_tx.setdefault(tx0, {}).setdefault((chrom_n, sgn), []).append((st, en))
            gene_by_tx[tx0] = gene0

    # Build CDSModel per transcript, and gene -> transcripts index
    tx_models: Dict[str, CDSModel] = {}
    gene2tx: Dict[str, List[str]] = {}

    for tx0, by_chrstrand in cds_by_tx.items():
        # choose chrom/strand with max CDS length
        best_key, best_len = None, -1
        for key, exs in by_chrstrand.items():
            L = sum(e - s + 1 for s, e in exs)
            if L > best_len:
                best_len, best_key = L, key
        if best_key is None:
            continue

        chrom_n, sgn = best_key
        exons = by_chrstrand[best_key]

        # transcription order
        exons = sorted(exons, key=lambda x: x[0], reverse=(sgn == -1))

        cum, running = [], 0
        for s, e in exons:
            cum.append(running)
            running += (e - s + 1)

        gene0 = gene_by_tx.get(tx0, "")
        tx_models[tx0] = CDSModel(
            chrom=chrom_n, strand=sgn, exons=exons, cumlen=cum, gene_id=gene0, tx_id=tx0
        )
        if gene0:
            gene2tx.setdefault(gene0, []).append(tx0)

    return tx_models, gene2tx

def pos_to_aapos(model: CDSModel, chrom: str, strand: int, pos1: int) -> Optional[int]:
    # pos1 is assumed 1-based genomic, consistent with GTF.
    # (We are NOT doing any pos+1 conversion trick here.)
    chrom_n = norm_chr(str(chrom))
    if chrom_n != model.chrom or int(strand) != model.strand:
        return None

    for (s, e), before in zip(model.exons, model.cumlen):
        if s <= pos1 <= e:
            if model.strand == 1:
                cds0 = before + (pos1 - s)
            else:
                cds0 = before + (e - pos1)
            return (cds0 // 3) + 1
    return None

def pick_first_id(val: str) -> str:
    if val is None:
        return ""
    s = str(val).strip()
    if not s or s == ".":
        return ""
    if ";" in s:
        for x in s.split(";"):
            x = x.strip()
            if x and x != ".":
                return x
        return ""
    return s

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet-in", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--parquet-out", required=True)
    ap.add_argument("--cache-pkl", default=None, help="Optional pickle cache for GTF models")
    args = ap.parse_args()

    parq_in = Path(args.parquet_in)
    gtf = Path(args.gtf)
    parq_out = Path(args.parquet_out)
    cache = Path(args.cache_pkl) if args.cache_pkl else None

    df = pd.read_parquet(parq_in)

    need = ["chr", "pos", "strand"]
    for c in need:
        if c not in df.columns:
            raise SystemExit(f"❌ Missing required column: {c}")

    tx_col = "transcript_id" if "transcript_id" in df.columns else None
    gene_col = "gene_id" if "gene_id" in df.columns else None
    if not tx_col and not gene_col:
        raise SystemExit("❌ Need at least transcript_id or gene_id in parquet to map to CDS")

    if cache and cache.exists():
        with open(cache, "rb") as f:
            tx_models, gene2tx = pickle.load(f)
    else:
        tx_models, gene2tx = build_models_from_gtf(gtf)
        if cache:
            cache.parent.mkdir(parents=True, exist_ok=True)
            with open(cache, "wb") as f:
                pickle.dump((tx_models, gene2tx), f)

    # Arrays
    chr_arr = df["chr"].astype(str).values
    pos_arr = pd.to_numeric(df["pos"], errors="coerce").fillna(-1).astype(int).values
    strand_arr = pd.to_numeric(df["strand"], errors="coerce").fillna(0).astype(int).values

    tx_arr = df[tx_col].astype(str).values if tx_col else None
    gene_arr = df[gene_col].astype(str).values if gene_col else None

    a1, a2, used_tx = [], [], []

    hits_tx = 0
    hits_gene = 0

    for i in range(len(df)):
        chrom = chr_arr[i]
        pos1 = pos_arr[i]          # assume this is already 1-based (no conversion)
        pos1_b = pos1 + 1          # second base of DBS
        sgn = int(strand_arr[i])

        if pos1 < 1:
            a1.append(None); a2.append(None); used_tx.append("")
            continue

        # --- try transcript_id first ---
        tx0 = strip_version(strip_version(pick_first_id(tx_arr[i])) if tx_arr is not None else "")
        model = tx_models.get(tx0) if tx0 else None

        aa1 = aa2 = None
        chosen = ""

        if model is not None:
            aa1 = pos_to_aapos(model, chrom, sgn, pos1)
            aa2 = pos_to_aapos(model, chrom, sgn, pos1_b)
            if aa1 is not None:
                hits_tx += 1
                chosen = model.tx_id

        # --- fallback: try all transcripts from gene_id ---
        if aa1 is None and gene_arr is not None:
            g0 = strip_version(pick_first_id(gene_arr[i]))
            for cand_tx in gene2tx.get(g0, []):
                m = tx_models.get(cand_tx)
                if m is None:
                    continue
                aa1t = pos_to_aapos(m, chrom, sgn, pos1)
                if aa1t is None:
                    continue
                aa2t = pos_to_aapos(m, chrom, sgn, pos1_b)
                aa1, aa2 = aa1t, aa2t
                hits_gene += 1
                chosen = m.tx_id
                break

        a1.append(aa1)
        a2.append(aa2)
        used_tx.append(chosen)

    df["aapos_primary"] = a1
    df["aapos_secondary"] = a2
    df["spans_two_codons"] = (df["aapos_primary"].notna() &
                              df["aapos_secondary"].notna() &
                              (df["aapos_primary"] != df["aapos_secondary"]))
    df["aapos_used_transcript"] = used_tx

    n = len(df)
    n1 = int(df["aapos_primary"].notna().sum())
    n2 = int(df["aapos_secondary"].notna().sum())
    nspan = int(df["spans_two_codons"].sum())

    print(f"✅ wrote {parq_out}")
    print(f"QC: aapos_primary present:   {n1}/{n} ({(n1/n*100):.1f}%)")
    print(f"QC: aapos_secondary present: {n2}/{n} ({(n2/n*100):.1f}%)")
    print(f"QC: spans_two_codons:        {nspan}/{n} ({(nspan/n*100):.1f}%)")
    print(f"QC: hits via transcript_id:  {hits_tx}/{n}")
    print(f"QC: hits via gene fallback:  {hits_gene}/{n}")
    print("Tip: if both hit counts are ~0, it’s almost certainly a coordinate convention mismatch OR chr naming mismatch.")

    parq_out.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(parq_out, index=False)

if __name__ == "__main__":
    main()
