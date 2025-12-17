#!/usr/bin/env python3
import argparse, re, pickle
from pathlib import Path
import pandas as pd

ATTR_RE = re.compile(r'(\S+)\s+"([^"]+)"')

def parse_attrs(s: str) -> dict:
    return {k:v for k,v in ATTR_RE.findall(s)}

def open_text(path: Path):
    if str(path).endswith(".gz"):
        import gzip
        return gzip.open(path, "rt")
    return open(path, "rt")

def load_cds_by_gene(gtf: Path, wanted_genes: set[str], cache: Path | None = None):
    """
    Returns:
      gene_id -> transcript_id -> {"strand": "+/-", "cds": [(start,end)...], "cds_len": int}
    CDS coords are 1-based inclusive, blocks stored in transcript order.
    """
    if cache and cache.exists():
        with open(cache, "rb") as f:
            return pickle.load(f)

    gene_map = {}
    with open_text(gtf) as fh:
        for line in fh:
            if not line or line[0] == "#":
                continue
            f = line.rstrip("\n").split("\t")
            if len(f) < 9:
                continue
            chrom, source, feature, start, end, score, strand, frame, attrs = f
            if feature != "CDS":
                continue
            a = parse_attrs(attrs)
            gid = a.get("gene_id")
            tid = a.get("transcript_id")
            if not gid or not tid or gid not in wanted_genes:
                continue
            s = int(start); e = int(end)
            gene_map.setdefault(gid, {}).setdefault(tid, {"strand": strand, "cds": [], "cds_len": 0})
            gene_map[gid][tid]["cds"].append((s, e))

    # sort blocks in transcript order + compute cds_len
    for gid, txs in gene_map.items():
        for tid, info in txs.items():
            blocks = info["cds"]
            if info["strand"] == "+":
                blocks.sort(key=lambda x: x[0])
            else:
                blocks.sort(key=lambda x: x[0], reverse=True)
            info["cds"] = blocks
            info["cds_len"] = sum(abs(e - s) + 1 for s, e in blocks)

    if cache:
        cache.parent.mkdir(parents=True, exist_ok=True)
        with open(cache, "wb") as f:
            pickle.dump(gene_map, f)

    return gene_map

def cds_offset_for_pos(pos1: int, blocks: list[tuple[int,int]], strand: str):
    offset = 0
    for (s, e) in blocks:
        lo, hi = (s, e) if s <= e else (e, s)
        if lo <= pos1 <= hi:
            if strand == "+":
                return offset + (pos1 - lo)
            else:
                return offset + (hi - pos1)
        offset += (hi - lo + 1)
    return None

def pick_transcript_containing_pos(txs: dict, pos1: int):
    """
    Choose among transcripts where pos is in CDS.
    Return (tid, info, off) or (None,None,None).
    Preference: transcript with largest CDS length.
    """
    candidates = []
    for tid, info in txs.items():
        off = cds_offset_for_pos(pos1, info["cds"], info["strand"])
        if off is not None:
            candidates.append((info["cds_len"], tid, info, off))
    if not candidates:
        return None, None, None
    candidates.sort(reverse=True)  # max cds_len first
    _, tid, info, off = candidates[0]
    return tid, info, off

def main():
    ap = argparse.ArgumentParser(description="Add aapos_primary/aapos_secondary using gene_id + GTF CDS.")
    ap.add_argument("--parquet-in", required=True)
    ap.add_argument("--gtf", required=True)
    ap.add_argument("--parquet-out", required=True)
    ap.add_argument("--cache-pkl", default=None, help="Optional cache for CDS intervals keyed by gene_id.")
    args = ap.parse_args()

    pin = Path(args.parquet_in)
    df = pd.read_parquet(pin)

    for col in ["gene_id", "pos", "chr", "strand"]:
        if col not in df.columns:
            raise SystemExit(f"❌ Missing required column '{col}' in {pin}")

    wanted_genes = set(df["gene_id"].dropna().astype(str).unique().tolist())
    gene_map = load_cds_by_gene(Path(args.gtf), wanted_genes, Path(args.cache_pkl) if args.cache_pkl else None)

    pos_arr = pd.to_numeric(df["pos"], errors="coerce").astype("Int64")

    cds_offset = []
    within = []
    a1 = []
    a2 = []
    chosen_tid = []

    for gid, pos1 in zip(df["gene_id"].astype(str), pos_arr):
        if pd.isna(pos1) or gid not in gene_map:
            cds_offset.append(None); within.append(None); a1.append(None); a2.append(None); chosen_tid.append(None)
            continue

        tid, info, off = pick_transcript_containing_pos(gene_map[gid], int(pos1))
        if off is None:
            cds_offset.append(None); within.append(None); a1.append(None); a2.append(None); chosen_tid.append(None)
            continue

        cds_offset.append(int(off))
        within_c = int(off % 3)
        within.append(within_c)
        aapos_primary = int(off // 3 + 1)
        a1.append(aapos_primary)
        chosen_tid.append(tid)

        # secondary base is pos+1
        off2 = cds_offset_for_pos(int(pos1) + 1, info["cds"], info["strand"])
        if off2 is None:
            a2.append(None)
        else:
            a2.append(int(off2 // 3 + 1))

    df["chosen_transcript_id"] = chosen_tid
    df["cds_offset"] = cds_offset
    df["within_codon"] = within
    df["aapos_primary"] = a1
    df["aapos_secondary"] = a2

    pout = Path(args.parquet_out)
    pout.parent.mkdir(parents=True, exist_ok=True)
    df.to_parquet(pout, index=False)

    n = len(df)
    n_ok = pd.Series(a1).notna().sum()
    n_sec = pd.Series(a2).notna().sum()
    n_span = (pd.Series(a2).notna() & (pd.Series(a2) != pd.Series(a1))).sum()

    print(f"✅ wrote {pout}")
    print(f"QC: aapos_primary present: {n_ok}/{n} ({n_ok/n*100:.1f}%)")
    print(f"QC: aapos_secondary present: {n_sec}/{n} ({n_sec/n*100:.1f}%)")
    print(f"QC: spans two codons: {n_span}/{n} ({n_span/n*100:.1f}%)")

if __name__ == "__main__":
    main()
