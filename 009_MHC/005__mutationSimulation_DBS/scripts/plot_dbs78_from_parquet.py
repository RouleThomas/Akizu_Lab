#!/usr/bin/env python3
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import argparse
from pathlib import Path

# Optional FASTA check
try:
    import pysam
except Exception:
    pysam = None

# ---- DBS78 (order = your COSMIC file) ----
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
CONTEXT2ID = {c:i for i,c in enumerate(CONTEXTS_78)}

GROUPS = ["AC","AT","CC","CG","CT","GC","TA","TC","TG","TT"]
GROUP_COLORS = {
    "AC": "#00BFC4","AT": "#1F77B4","CC": "#8BC34A","CG": "#2E7D32","CT": "#F48FB1",
    "GC": "#E53935","TA": "#FFB74D","TC": "#FB8C00","TG": "#B39DDB","TT": "#6A1B9A",
}

def revcomp(seq: str) -> str:
    return seq.translate(str.maketrans("ACGT","TGCA"))[::-1]

def main(parquet_path, output_pdf, sample_name, fasta_path=None):
    df = pd.read_parquet(parquet_path)

    if "context_id" not in df.columns:
        raise ValueError("Parquet must contain a 'context_id' column (0–77).")

    # ---------- Mode A: directly from context_id ----------
    ctx_from_id = [CONTEXTS_78[int(i)] for i in df["context_id"].tolist()]
    counts = pd.Series(ctx_from_id).value_counts().reindex(CONTEXTS_78, fill_value=0)

    # ---------- (Optional) Mode B: FASTA sanity check ----------
    mismatches = 0
    if fasta_path:
        if pysam is None:
            print("pysam not available; skipping FASTA checks.")
        else:
            fasta = pysam.FastaFile(fasta_path)
            for chrom, pos, cid in zip(df["chr"], df["pos"], df["context_id"]):
                cid = int(cid)
                ref_can, alt_can = CONTEXTS_78[cid].split(">")
                try:
                    ref_gen = fasta.fetch(str(chrom), int(pos)-1, int(pos)+1).upper()
                except Exception:
                    mismatches += 1
                    continue
                # accept either orientation (simulation may store canonical only)
                if ref_gen != ref_can and ref_gen != revcomp(ref_can):
                    mismatches += 1
            print(f"FASTAsanity: {mismatches} / {len(df)} sites not matching canonical or RC ref dinucleotide.")

    # ---------- Plot ----------
    percent = 100 * counts / max(counts.sum(), 1)
    bar_colors = []
    groups_idx = {g: [] for g in GROUPS}
    for i, ctx in enumerate(CONTEXTS_78):
        ref2 = ctx.split(">")[0]
        grp = ref2[:2]
        bar_colors.append(GROUP_COLORS.get(grp, "#999999"))
        groups_idx[grp].append(i)

    fig, ax = plt.subplots(figsize=(22, 6))
    ax.bar(range(78), percent.values, color=bar_colors, width=0.95)
    ax.set_xticks(range(78))
    ax.set_xticklabels(CONTEXTS_78, rotation=90, fontsize=6)
    ax.set_ylabel("Percentage of Double Base Substitutions")
    ax.set_title(sample_name)

    y_top = float(percent.max()) if len(percent) else 0.0
    ax.set_ylim(0, max(y_top * 1.15, 1))
    band_y = y_top * 1.05
    for grp, idxs in groups_idx.items():
        x0, x1 = min(idxs)-0.5, max(idxs)+0.5
        ax.axvspan(x0, x1, color=GROUP_COLORS[grp], alpha=0.10)
        ax.text((x0+x1)/2, band_y, f"{grp}>NN", ha="center", va="bottom",
                fontsize=10, color=GROUP_COLORS[grp], fontweight="bold")

    fig.tight_layout()
    fig.savefig(output_pdf, dpi=300)
    print(f"📊 DBS78 plot saved: {output_pdf}")

if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--parquet", required=True)
    ap.add_argument("--output-pdf", required=True)
    ap.add_argument("--sample-name", required=True)
    ap.add_argument("--fasta", help="Optional: GRCh38 FASTA for sanity checks")
    args = ap.parse_args()
    main(args.parquet, args.output_pdf, args.sample_name, args.fasta)
