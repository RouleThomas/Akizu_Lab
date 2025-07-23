#!/bin/bash
#SBATCH --job-name=exome_tab
#SBATCH --array=1-22
#SBATCH --mem=50G
#SBATCH --time=2:00:00
#SBATCH --output=logs/sim_%A_%a.out
#SBATCH --error=logs/sim_%A_%a.err


# Chromosome to process (1–22 from SLURM array)
CHR=${SLURM_ARRAY_TASK_ID}

# Run Python script to generate per-nucleotide .parquet for one chromosome
python scripts/build_exome_table.py \
  --cds-fa gtf/cds.fa \
  --chrom chr${CHR} \
  --output parquet/chr${CHR}.parquet


