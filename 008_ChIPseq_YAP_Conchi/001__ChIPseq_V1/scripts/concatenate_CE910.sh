#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00

cat input_raw/CE9_S9_L001_R1_001.fastq.gz input_raw/CE9_S9_L002_R1_001.fastq.gz input_raw/CE9_S9_L003_R1_001.fastq.gz input_raw/CE9_S9_L004_R1_001.fastq.gz > input_raw/CE9_S9_L14_R1_001.fastq.gz
cat input_raw/CE9_S9_L001_R2_001.fastq.gz input_raw/CE9_S9_L002_R2_001.fastq.gz input_raw/CE9_S9_L003_R2_001.fastq.gz input_raw/CE9_S9_L004_R2_001.fastq.gz > input_raw/CE9_S9_L14_R2_001.fastq.gz

cat input_raw/CE10_S9_L001_R1_001.fastq.gz input_raw/CE10_S9_L002_R1_001.fastq.gz input_raw/CE10_S9_L003_R1_001.fastq.gz input_raw/CE10_S9_L004_R1_001.fastq.gz > input_raw/CE10_S9_L14_R1_001.fastq.gz
cat input_raw/CE10_S9_L001_R2_001.fastq.gz input_raw/CE10_S9_L002_R2_001.fastq.gz input_raw/CE10_S9_L003_R2_001.fastq.gz input_raw/CE10_S9_L004_R2_001.fastq.gz > input_raw/CE10_S9_L14_R2_001.fastq.gz


