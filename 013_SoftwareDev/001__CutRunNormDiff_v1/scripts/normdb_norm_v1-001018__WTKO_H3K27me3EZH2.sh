#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=5



set -euo pipefail  # if any bug, stop the script



python -u scripts/normdb_norm_v1.py normalize \
  --meta meta/samples-001018__WTKO_H3K27me3EZH2.tsv \
  --outdir output/normdb_norm-001018__WTKO_H3K27me3EZH2 \
  --blacklist meta/hg38-blacklist.v2.bed \
  --chrom-sizes meta/GRCh38_chrom_sizes.tab \
  --threads 5 \
  --mode PE \
  --reference auto






