#!/usr/bin/env python3
import argparse
from pathlib import Path
import pandas as pd
import numpy as np
import pysam

# DBS78 contexts (order must match your simulation)
CONTEXTS_78 = [
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
def revcomp(s: str) -> str:
    return s.translate(str.maketrans("ACGT","TGCA"))[::-1]

def main(parquet_in, fasta_path, parquet_out, tsv_out=None, vcf_out=None, preview_n=1000):
    df = pd.read_parquet(parquet_in)
    required = {"chr","pos","context_id"}
    if not required.issubset(df.columns):
        missing = ", ".join(sorted(required - set(df.columns)))
        raise SystemExit(f"Parquet missing required columns: {missing}")

    # normalize types
    df["context_id"] = df["context_id"].astype(int)
    if "strand" in df.columns:
        # convert strand to +/-
        strand_sign = df["strand"].astype(int).map({1:"+", -1:"-"}).fillna(".")
    else:
        strand_sign = pd.Series(["."]*len(df), index=df.index)

    fasta = pysam.FastaFile(str(fasta_path))

    ref2_list, alt2_list, ok_list, note_list = [], [], [], []
    ctx_list, orient_list, txref_list, txalt_list = [], [], [], []

    for i, (chrom, pos, cid, sgn) in enumerate(zip(df["chr"], df["pos"], df["context_id"], strand_sign)):
        pos = int(pos); cid = int(cid)
        ctx_can = CONTEXTS_78[cid]
        ref_can, alt_can = ctx_can.split(">")

        # fetch genome dinucleotide at leftmost base (POS)
        try:
            ref_gen = fasta.fetch(str(chrom), pos-1, pos+1).upper()
        except Exception as e:
            ref2_list.append(None); alt2_list.append(None)
            ok_list.append(False); note_list.append(f"fetch_fail:{e}")
            ctx_list.append(ctx_can); orient_list.append("NA"); txref_list.append(None); txalt_list.append(None)
            continue

        # decide orientation and set genome-oriented REF/ALT
        if ref_gen == ref_can:
            orient = "canon"; ref2, alt2 = ref_can, alt_can
        elif ref_gen == revcomp(ref_can):
            orient = "rc";    ref2, alt2 = ref_gen, revcomp(alt_can)
        else:
            ref2, alt2 = None, None
            ok_list.append(False); note_list.append(f"ref_mismatch:genome={ref_gen},canon={ref_can}")
            ctx_list.append(ctx_can); orient_list.append("NA"); txref_list.append(None); txalt_list.append(None)
            continue

        # transcript-oriented alleles (if strand known)
        if sgn in {"+","-"}:
            txref = ref2 if sgn == "+" else revcomp(ref2)
            txalt = alt2 if sgn == "+" else revcomp(alt2)
        else:
            txref = txalt = None

        ref2_list.append(ref2); alt2_list.append(alt2)
        ok_list.append(True); note_list.append("canonical" if orient=="canon" else "reverse_complement")
        ctx_list.append(ctx_can); orient_list.append(orient)
        txref_list.append(txref); txalt_list.append(txalt)

    # attach to dataframe
    df["ref2"] = ref2_list
    df["alt2"] = alt2_list
    df["match_ok"] = ok_list
    df["note"] = note_list
    df["CTX"] = ctx_list
    df["ORIENT"] = orient_list
    df["STRAND_SIGN"] = strand_sign
    df["TXREF"] = txref_list
    df["TXALT"] = txalt_list

    # save updated parquet
    df.to_parquet(parquet_out, compression="zstd")
    print(f"✅ wrote Parquet: {parquet_out}  rows={len(df)}  ok={int(df.match_ok.sum())}  bad={int((~df.match_ok).sum())}")

    # preview TSV (includes STRAND + context/orientation + TXREF/TXALT)
    if tsv_out:
        preview_cols = ["chr","pos","ref2","alt2","context_id","CTX","ORIENT","STRAND_SIGN","TXREF","TXALT",
                        "signature","rep","seed","match_ok","note","gene_id","transcript_id"]
        preview_cols = [c for c in preview_cols if c in df.columns]
        df.head(preview_n)[preview_cols].to_csv(tsv_out, sep="\t", index=False)
        print(f"📝 wrote preview TSV: {tsv_out} (first {min(preview_n, len(df))} rows)")

    # minimal VCF (only ok rows)
    if vcf_out:
        with open(vcf_out, "w") as f:
            f.write("##fileformat=VCFv4.2\n")
            f.write("##INFO=<ID=CTX,Number=1,Type=String,Description=\"Canonical pyrimidine-left DBS context (ref>alt)\">\n")
            f.write("##INFO=<ID=ORIENT,Number=1,Type=String,Description=\"Allele orientation relative to canonical (canon|rc)\">\n")
            f.write("##INFO=<ID=STRAND,Number=1,Type=String,Description=\"Transcript strand if present (+/-/. )\">\n")
            f.write("##INFO=<ID=TXREF,Number=1,Type=String,Description=\"Allele on transcript strand (if strand known)\">\n")
            f.write("##INFO=<ID=TXALT,Number=1,Type=String,Description=\"Allele on transcript strand (if strand known)\">\n")
            f.write("#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n")
            for chrom, pos, ok, ref2, alt2, ctx, orient, sgn, txr, txa in \
                zip(df["chr"], df["pos"], df["match_ok"], df["ref2"], df["alt2"], df["CTX"], df["ORIENT"], df["STRAND_SIGN"], df["TXREF"], df["TXALT"]):
                if not ok or ref2 is None or alt2 is None:
                    continue
                info = [f"CTX={ctx}", f"ORIENT={orient}", f"STRAND={sgn}"]
                if txr is not None and txa is not None:
                    info += [f"TXREF={txr}", f"TXALT={txa}"]
                f.write(f"{chrom}\t{int(pos)}\t.\t{ref2}\t{alt2}\t.\t.\t{';'.join(info)}\n")
        print(f"🧬 wrote VCF: {vcf_out}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser(description="Attach genome-oriented REF/ALT to DBS parquet; emit preview TSV and VCF for IGV.")
    ap.add_argument("--parquet-in", required=True)
    ap.add_argument("--fasta", required=True, help="Reference FASTA (with .fai)")
    ap.add_argument("--parquet-out", required=True)
    ap.add_argument("--tsv-out", default=None)
    ap.add_argument("--vcf-out", default=None)
    ap.add_argument("--preview", type=int, default=1000)
    args = ap.parse_args()
    main(Path(args.parquet_in), Path(args.fasta), Path(args.parquet_out),
         Path(args.tsv_out) if args.tsv_out else None,
         Path(args.vcf_out) if args.vcf_out else None,
         preview_n=args.preview)
