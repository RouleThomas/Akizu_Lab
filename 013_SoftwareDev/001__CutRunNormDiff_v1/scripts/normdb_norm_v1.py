#!/usr/bin/env python3
import argparse
import os
import sys
import shutil
import subprocess
import time
from datetime import datetime

import pandas as pd
import numpy as np

REQUIRED_COLS = ["sample_id", "bam", "condition", "target"]


# -----------------------------
# Helpers
# -----------------------------
def now():
    return datetime.now().strftime("%H:%M:%S")


def log_step(msg, start_time=None):
    if start_time:
        elapsed = time.time() - start_time
        print(f"[{now()}] {msg} (done in {elapsed:.1f}s)")
    else:
        print(f"[{now()}] {msg}")


def run_cmd(cmd, log_fh):
    cmd_str = " ".join(cmd)
    log_fh.write(f"\n[{datetime.now().isoformat(timespec='seconds')}] CMD: {cmd_str}\n")
    log_fh.flush()

    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    log_fh.write(res.stdout)
    log_fh.flush()

    if res.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd_str}\n{res.stdout}")


def ensure_exe(name):
    if shutil.which(name) is None:
        raise RuntimeError(f"Required executable not found in PATH: {name}")


# -----------------------------
# Core logic
# -----------------------------
def read_meta(meta_path):
    df = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        raise ValueError(f"Metadata missing required columns: {missing}")

    if df["sample_id"].duplicated().any():
        raise ValueError("Duplicate sample_id detected")

    for bam in df["bam"]:
        if not bam.endswith(".bam") or not os.path.exists(bam):
            raise ValueError(f"Invalid BAM file: {bam}")

    return df


def find_local_maxima_bedgraph(bedgraph_path, out_bed_path):
    data = pd.read_csv(
        bedgraph_path, sep="\t", header=None,
        names=["chrom", "start", "end", "score"]
    )

    scores = data["score"].to_numpy()
    if len(scores) < 3:
        open(out_bed_path, "w").close()
        return 0

    idx = []
    for i in range(1, len(scores) - 1):
        if scores[i] > scores[i - 1] and scores[i] > scores[i + 1]:
            idx.append(i)

    data.iloc[idx].to_csv(out_bed_path, sep="\t", index=False, header=False)
    return len(idx)


def percentile_99_from_maxima(path):
    df = pd.read_csv(path, sep="\t", header=None, names=["c", "s", "e", "v"])
    if df.empty:
        return None
    return float(np.percentile(df["v"], 99))


# -----------------------------
# Main
# -----------------------------
def main():
    t0 = time.time()

    parser = argparse.ArgumentParser("normdb")
    sub = parser.add_subparsers(dest="command", required=True)

    p = sub.add_parser("normalize")
    p.add_argument("--meta", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--blacklist")
    p.add_argument("--chrom-sizes", required=True)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--mode", choices=["PE", "SE"], required=True)
    p.add_argument("--se-fragment-length", type=int, default=200)
    p.add_argument("--reference", default="auto")

    args = parser.parse_args()

    # checks
    for exe in ["bamCoverage", "bigWigToBedGraph", "bedGraphToBigWig", "bedtools"]:
        ensure_exe(exe)

    os.makedirs(args.outdir, exist_ok=True)

    # output structure
    raw_bw = f"{args.outdir}/01_raw_bigwig"
    bg = f"{args.outdir}/02_bedgraph"
    bg_bl = f"{args.outdir}/03_bedgraph_blacklist"
    maxima = f"{args.outdir}/04_local_maxima"
    norm_bg = f"{args.outdir}/05_normalized_bedgraph"
    norm_bw = f"{args.outdir}/06_normalized_bigwig"

    for d in [raw_bw, bg, bg_bl, maxima, norm_bg, norm_bw]:
        os.makedirs(d, exist_ok=True)

    log_path = f"{args.outdir}/normdb_normalize.log"

    with open(log_path, "w") as log:
        log.write(f"Started: {datetime.now()}\n")
        log.write(f"Args: {vars(args)}\n")

        # -----------------------------
        # Load metadata
        # -----------------------------
        step = time.time()
        log_step("Reading metadata")
        meta = read_meta(args.meta)
        log_step("Metadata loaded", step)

        # -----------------------------
        # BAM → bigWig → bedGraph → blacklist
        # -----------------------------
        step = time.time()
        log_step("Converting BAM → bigWig → bedGraph")
        for _, r in meta.iterrows():
            sid = r.sample_id

            bw = f"{raw_bw}/{sid}.bw"
            bgf = f"{bg}/{sid}.bedGraph"
            bgf_bl = f"{bg_bl}/{sid}.bedGraph"

            cmd = [
                "bamCoverage",
                "--bam", r.bam,
                "--outFileName", bw,
                "--outFileFormat", "bigwig",
                "--binSize", "1",
                "--numberOfProcessors", str(args.threads),
                "--scaleFactor", "1",
            ]
            if args.mode == "PE":
                cmd += ["--extendReads"]
            else:
                cmd += ["--extendReads", str(args.se_fragment_length)]

            run_cmd(cmd, log)
            run_cmd(["bigWigToBedGraph", bw, bgf], log)

            if args.blacklist:
                with open(bgf_bl, "w") as out:
                    subprocess.run(
                        ["bedtools", "intersect", "-v", "-a", bgf, "-b", args.blacklist],
                        stdout=out, check=True
                    )
            else:
                shutil.copy(bgf, bgf_bl)

        log_step("Raw signal generation finished", step)

        # -----------------------------
        # Local maxima
        # -----------------------------
        step = time.time()
        log_step("Identifying local maxima")
        for sid in meta.sample_id:
            find_local_maxima_bedgraph(
                f"{bg_bl}/{sid}.bedGraph",
                f"{maxima}/{sid}.local_maxima.bed"
            )
        log_step("Local maxima identified", step)

        # -----------------------------
        # P99 + scaling
        # -----------------------------
        step = time.time()
        log_step("Computing 99th percentiles and scaling factors")
        p99 = {}
        for sid in meta.sample_id:
            val = percentile_99_from_maxima(f"{maxima}/{sid}.local_maxima.bed")
            if val is None:
                raise RuntimeError(f"No maxima for {sid}")
            p99[sid] = val

        ref_by_target = {}
        for t in meta.target.unique():
            ref = meta[meta.target == t].iloc[0].sample_id
            ref_by_target[t] = ref

        scaling = {}
        for _, r in meta.iterrows():
            scaling[r.sample_id] = p99[ref_by_target[r.target]] / p99[r.sample_id]

        log_step("Scaling factors computed", step)

        # -----------------------------
        # Normalize + bigWig
        # -----------------------------
        step = time.time()
        log_step("Generating normalized bigWig files")
        for sid in meta.sample_id:
            df = pd.read_csv(f"{bg_bl}/{sid}.bedGraph", sep="\t", header=None)
            df[3] *= scaling[sid]

            norm_bgf = f"{norm_bg}/{sid}.norm99.bedGraph"
            norm_bgf_sorted = f"{norm_bg}/{sid}.norm99.sorted.bedGraph"
            norm_bwf = f"{norm_bw}/{sid}.norm99.bw"

            df.to_csv(norm_bgf, sep="\t", index=False, header=False)
            subprocess.run(
                ["bedtools", "sort", "-i", norm_bgf],
                stdout=open(norm_bgf_sorted, "w"),
                check=True
            )
            run_cmd(["bedGraphToBigWig", norm_bgf_sorted, args.chrom_sizes, norm_bwf], log)

        log_step("Normalization finished", step)

        log.write(f"Finished: {datetime.now()}\n")

    print(f"\n✔ NORMDB normalize finished in {(time.time() - t0)/60:.1f} minutes")
    print(f"Output directory: {args.outdir}")
    print(f"Log file: {log_path}")


if __name__ == "__main__":
    main()
