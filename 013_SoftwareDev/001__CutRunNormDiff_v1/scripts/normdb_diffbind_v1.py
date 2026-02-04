#!/usr/bin/env python3
import argparse
import os
import sys
import shutil
import subprocess
import time
from datetime import datetime

import pandas as pd


REQUIRED_COLS = ["sample_id", "bam", "condition", "target"]


# -----------------------------
# Pretty logging + timing
# -----------------------------
def now():
    return datetime.now().strftime("%H:%M:%S")


def log(msg):
    print(f"[{now()}] {msg}", flush=True)


def log_done(msg, tstart):
    log(f"{msg} (done in {time.time() - tstart:.1f}s)")


def die(msg):
    print(f"ERROR: {msg}", file=sys.stderr)
    sys.exit(1)


def ensure_exe(name):
    if shutil.which(name) is None:
        die(f"Required executable not found in PATH: {name}")


def run_cmd(cmd, log_fh=None):
    cmd_str = " ".join(cmd)
    if log_fh:
        log_fh.write(f"\n[{datetime.now().isoformat(timespec='seconds')}] CMD: {cmd_str}\n")
        log_fh.flush()

    res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    if log_fh:
        log_fh.write(res.stdout)
        log_fh.flush()

    if res.returncode != 0:
        raise RuntimeError(f"Command failed:\n{cmd_str}\n{res.stdout}")
    return res.stdout


# -----------------------------
# IO helpers
# -----------------------------
def read_meta(meta_path):
    if not os.path.exists(meta_path):
        die(f"Metadata file not found: {meta_path}")

    df = pd.read_csv(meta_path, sep="\t", dtype=str).fillna("")
    missing = [c for c in REQUIRED_COLS if c not in df.columns]
    if missing:
        die(f"Metadata missing required columns: {missing} (required: {REQUIRED_COLS})")

    if (df["sample_id"].str.strip() == "").any():
        die("Metadata contains empty sample_id.")
    if df["sample_id"].duplicated().any():
        dups = df[df["sample_id"].duplicated()]["sample_id"].tolist()
        die(f"Metadata contains duplicate sample_id(s): {dups}")

    return df


def parse_contrast(s):
    # expected: condition:KO:WT
    parts = s.split(":")
    if len(parts) != 3:
        die("--contrast must look like: condition:KO:WT")
    return parts[0], parts[1], parts[2]


def read_computeMatrix_vector(txt_path):
    """
    Your original read.delim(..., skip=3) suggested the computeMatrix .txt has 3 header lines.
    With regionBodyLength=100 and binSize=100 => exactly ONE bin per region.
    After skipping headers, we take the first column as the summed signal per region.
    """
    if not os.path.exists(txt_path):
        die(f"computeMatrix matrix file not found: {txt_path}")

    # robust parsing: skip comment lines + potential header
    # We mimic: read.delim(..., header=FALSE, skip=3)
    df = pd.read_csv(txt_path, sep="\t", header=None, skiprows=3)
    if df.shape[1] < 1:
        die(f"computeMatrix matrix file seems empty/unexpected: {txt_path}")

    vec = df.iloc[:, 0].astype(float).tolist()
    return vec


def bed_to_region_ids(bed_path):
    """
    computeMatrix --outFileSortedRegions produces a BED-like file with a header in some versions.
    We'll parse the first 3 columns, skipping non-data lines.
    """
    if not os.path.exists(bed_path):
        die(f"computeMatrix sorted regions file not found: {bed_path}")

    rows = []
    with open(bed_path) as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            # Some outputs have column names like "X.chrom" in first line; detect non-bed
            if parts[0].lower() in ("chrom", "chr", "x.chrom", "x.chromosome"):
                continue
            if len(parts) < 3:
                continue
            chrom = parts[0]
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            rows.append((chrom, start, end))

    if not rows:
        die(f"Could not parse any regions from: {bed_path}")

    return [f"{c}_{s}_{e}" for (c, s, e) in rows]


# -----------------------------
# Main
# -----------------------------
def main():
    t0 = time.time()

    ap = argparse.ArgumentParser(prog="normdb diffbind", description="Diff binding on user BED regions (computeMatrix-based)")
    ap.add_argument("--meta", required=True, help="samples.tsv with columns: sample_id bam condition target")
    ap.add_argument("--regions", required=True, help="BED file of regions (chr start end)")
    ap.add_argument("--bigwig-dir", required=True, help="Directory containing normalized bigWigs")
    ap.add_argument("--contrast", required=True, help="e.g. condition:KO:WT")
    ap.add_argument("--outdir", required=True, help="Output directory")
    ap.add_argument("--alpha", type=float, default=0.05, help="Used for plot highlighting only")
    ap.add_argument("--lfc", type=float, default=0.0, help="Used for plot highlighting only (log2FC threshold)")
    ap.add_argument("--target", default=None, help="If meta contains multiple targets, specify which one (e.g. H3K27me3)")

    args = ap.parse_args()

    # Required executables
    ensure_exe("computeMatrix")
    ensure_exe("Rscript")

    if not os.path.exists(args.regions):
        die(f"--regions not found: {args.regions}")
    if not os.path.isdir(args.bigwig_dir):
        die(f"--bigwig-dir not found or not a directory: {args.bigwig_dir}")

    os.makedirs(args.outdir, exist_ok=True)

    # Subfolders (like your old layout)
    cm_dir = os.path.join(args.outdir, "01_computeMatrix")
    os.makedirs(cm_dir, exist_ok=True)

    log_path = os.path.join(args.outdir, "normdb_diffbind.log")
    with open(log_path, "w") as log_fh:
        log_fh.write(f"Started: {datetime.now().isoformat(timespec='seconds')}\n")
        log_fh.write(f"Args: {vars(args)}\n")

        # -------------------------
        # Load + filter meta
        # -------------------------
        step = time.time()
        log("Reading metadata")
        meta = read_meta(args.meta)

        targets = sorted(meta["target"].unique().tolist())
        if args.target is None:
            if len(targets) != 1:
                die(f"Meta contains multiple targets {targets}. Please provide --target <one>.")
            target = targets[0]
        else:
            if args.target not in targets:
                die(f"--target {args.target} not found in meta targets: {targets}")
            target = args.target

        meta_t = meta.loc[meta["target"] == target].copy()
        if meta_t.shape[0] < 2:
            die(f"Not enough samples for target={target} (need >=2).")

        contrast_var, contrast_test, contrast_ref = parse_contrast(args.contrast)
        if contrast_var not in meta_t.columns:
            die(f"Contrast variable '{contrast_var}' not found in meta columns: {list(meta_t.columns)}")

        levels = set(meta_t[contrast_var].tolist())
        if contrast_test not in levels or contrast_ref not in levels:
            die(f"Meta does not contain both contrast levels: {contrast_test}, {contrast_ref}. Found: {sorted(levels)}")

        log(f"Using target={target} with n={meta_t.shape[0]} samples; contrast={contrast_var}:{contrast_test}:{contrast_ref}")
        log_done("Metadata loaded", step)

        # -------------------------
        # Collect bigWigs
        # Expect: <bigwig-dir>/<sample_id>.norm99.bw
        # -------------------------
        step = time.time()
        log("Collecting bigWig files")
        bw_paths = []
        bw_labels = []
        missing = []

        # Preserve metadata order
        for _, r in meta_t.iterrows():
            sid = r["sample_id"]
            bw = os.path.join(args.bigwig_dir, f"{sid}.norm99.bw")
            if not os.path.exists(bw):
                missing.append(bw)
            else:
                bw_paths.append(bw)
                bw_labels.append(sid)

        if missing:
            die("Missing bigWig(s):\n" + "\n".join(missing))

        log_done(f"Found {len(bw_paths)} bigWigs", step)

        # -------------------------
        # Step 1: computeMatrix per sample (EXACT settings)
        # -------------------------
        step = time.time()
        log("computeMatrix: counting signal in regions (exact scale-regions method)")

        matrix_txts = {}
        sorted_beds = {}
        for sid, bw in zip(bw_labels, bw_paths):
            # Match your naming vibe (short, consistent)
            base = f"{sid}-{target}"
            out_gz = os.path.join(cm_dir, f"{base}.gz")
            out_txt = os.path.join(cm_dir, f"{base}.txt")
            out_bed = os.path.join(cm_dir, f"{base}.bed")

            cmd = [
                "computeMatrix", "scale-regions",
                "-S", bw,
                "-R", args.regions,
                "--outFileName", out_gz,
                "--outFileNameMatrix", out_txt,
                "--outFileSortedRegions", out_bed,
                "--missingDataAsZero",
                "--averageTypeBins", "sum",
                "--binSize", "100",
                "--regionBodyLength", "100"
            ]
            run_cmd(cmd, log_fh=log_fh)

            matrix_txts[sid] = out_txt
            sorted_beds[sid] = out_bed

        log_done("computeMatrix finished for all samples", step)

        # -------------------------
        # Step 2: build count matrix (one value per region)
        # Use region order from FIRST sample’s outFileSortedRegions
        # -------------------------
        step = time.time()
        log("Building region count matrix from computeMatrix outputs")

        first_sid = bw_labels[0]
        region_ids = bed_to_region_ids(sorted_beds[first_sid])
        n_regions = len(region_ids)

        counts = pd.DataFrame({"region_id": region_ids})

        for sid in bw_labels:
            vec = read_computeMatrix_vector(matrix_txts[sid])
            if len(vec) != n_regions:
                die(
                    f"Region count mismatch for {sid}: "
                    f"{len(vec)} values in matrix but {n_regions} regions in sorted BED. "
                    f"This should not happen; check regions.bed and computeMatrix output."
                )
            # DESeq2 expects integer-ish counts: round
            counts[sid] = pd.Series(vec).round().astype(int)

        counts_path = os.path.join(args.outdir, "counts_matrix.tsv")
        counts.to_csv(counts_path, sep="\t", index=False)

        # coldata.tsv
        meta_t = meta_t.set_index("sample_id").loc[bw_labels].reset_index()
        coldata_path = os.path.join(args.outdir, "coldata.tsv")
        meta_t[["sample_id", contrast_var]].to_csv(coldata_path, sep="\t", index=False)

        log_done("Count matrix + coldata written", step)

        # -------------------------
        # Step 3: DESeq2 + plots
        # -------------------------
        step = time.time()
        log("Running DESeq2 + plots (volcano, MA, sample correlation heatmap)")

        results_path = os.path.join(args.outdir, "results_all_regions.tsv")
        volcano_pdf = os.path.join(args.outdir, "volcano.pdf")
        ma_pdf = os.path.join(args.outdir, "MA.pdf")
        corr_pdf = os.path.join(args.outdir, "sample_correlation_heatmap.pdf")
        r_script_path = os.path.join(args.outdir, "run_deseq2.R")

        r_script = f"""
suppressPackageStartupMessages({{
  library(DESeq2)
  library(tidyverse)
  library(EnhancedVolcano)
}})

alpha <- {args.alpha}
lfc_thr <- {args.lfc}

counts_file <- "{counts_path}"
coldata_file <- "{coldata_path}"

out_results <- "{results_path}"
out_volcano <- "{volcano_pdf}"
out_ma <- "{ma_pdf}"
out_corr <- "{corr_pdf}"

contrast_var <- "{contrast_var}"
test_level <- "{contrast_test}"
ref_level <- "{contrast_ref}"

counts_df <- read.delim(counts_file, check.names=FALSE)
stopifnot("region_id" %in% colnames(counts_df))

count_mat <- counts_df %>%
  select(-region_id) %>%
  as.matrix()

rownames(count_mat) <- counts_df$region_id
mode(count_mat) <- "integer"

coldata <- read.delim(coldata_file, check.names=FALSE)
stopifnot("sample_id" %in% colnames(coldata))
stopifnot(contrast_var %in% colnames(coldata))

# Ensure sample order matches matrix columns
stopifnot(all(coldata$sample_id == colnames(count_mat)))

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = as.data.frame(coldata %>% column_to_rownames("sample_id")),
  design = as.formula(paste0("~ ", contrast_var))
)

dds[[contrast_var]] <- relevel(dds[[contrast_var]], ref = ref_level)
dds <- DESeq(dds)

coef_name <- resultsNames(dds)[2]
use_apeglm <- requireNamespace("apeglm", quietly=TRUE)

if (use_apeglm) {{
  res <- lfcShrink(dds, coef=coef_name, type="apeglm")
}} else {{
  res <- lfcShrink(dds, coef=coef_name, type="normal")
}}

res_df <- as.data.frame(res) %>%
  rownames_to_column("region_id") %>%
  mutate(
    chr = sub("_(.*)$", "", region_id),
    start = as.integer(sub("^.*?_(\\d+)_.*$", "\\\\1", region_id)),
    end = as.integer(sub("^.*?_\\d+_(\\d+)$", "\\\\1", region_id))
  ) %>%
  relocate(chr, start, end, region_id)

# Write ALL regions
write.table(res_df, file=out_results, sep="\\t", quote=FALSE, row.names=FALSE)

# Volcano (highlighting only)
res_plot <- res_df %>%
  mutate(
    padj_plot = ifelse(is.na(padj), 1, padj),
    sig = (padj_plot < alpha) & (abs(log2FoldChange) >= lfc_thr),
    direction = case_when(
      sig & log2FoldChange > 0 ~ "Gain",
      sig & log2FoldChange < 0 ~ "Loss",
      TRUE ~ "NS"
    )
  )

keyvals <- ifelse(res_plot$direction == "Gain", "orange",
           ifelse(res_plot$direction == "Loss", "skyblue", "grey70"))
names(keyvals)[keyvals == "orange"] <- "Gain"
names(keyvals)[keyvals == "skyblue"] <- "Loss"
names(keyvals)[keyvals == "grey70"] <- "NS"

pdf(out_volcano, width=4.5, height=4.5)
EnhancedVolcano(
  res_plot,
  lab = rep("", nrow(res_plot)),
  x = "log2FoldChange",
  y = "padj_plot",
  pCutoff = alpha,
  FCcutoff = lfc_thr,
  colCustom = keyvals,
  pointSize = 1.0,
  labSize = 2.0,
  title = paste0(test_level, " vs ", ref_level, " (", contrast_var, ")"),
  subtitle = "computeMatrix scale-regions (sum, 1 bin/region) + DESeq2",
  legendPosition = "right"
) + theme_bw()
dev.off()

# MA plot
pdf(out_ma, width=4.5, height=4.0)
plotMA(res, alpha=alpha, main=paste0("MA: ", test_level, " vs ", ref_level))
abline(h=c(-lfc_thr, lfc_thr), lty=2)
dev.off()

# Sample correlation heatmap (vst)
vsd <- vst(dds, blind=FALSE)
mat <- assay(vsd)
cormat <- cor(mat, method="pearson")

pdf(out_corr, width=5.5, height=5.5)
heatmap(
  cormat,
  symm=TRUE,
  margins=c(8,8),
  main="Sample correlation (vst counts)"
)
dev.off()
"""
        with open(r_script_path, "w") as f:
            f.write(r_script)

        run_cmd(["Rscript", r_script_path], log_fh=log_fh)

        log_done("DESeq2 + plots finished", step)

    log(f"All done. Total runtime: {(time.time()-t0)/60:.1f} minutes")
    print(f"Outputs in: {args.outdir}")
    print(f"Log: {log_path}")


if __name__ == "__main__":
    main()
