#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00


cat input_raw/CE7_S11_L001_R1_001.fastq.gz input_raw/CE7_S11_L002_R1_001.fastq.gz input_raw/CE7_S11_L003_R1_001.fastq.gz input_raw/CE7_S11_L004_R1_001.fastq.gz > input_raw/CE7_S11_L14_R1_001.fastq.gz
cat input_raw/CE7_S11_L001_R2_001.fastq.gz input_raw/CE7_S11_L002_R2_001.fastq.gz input_raw/CE7_S11_L003_R2_001.fastq.gz input_raw/CE7_S11_L004_R2_001.fastq.gz > input_raw/CE7_S11_L14_R2_001.fastq.gz

cat input_raw/CE8_S8_L001_R1_001.fastq.gz input_raw/CE8_S8_L002_R1_001.fastq.gz input_raw/CE8_S8_L003_R1_001.fastq.gz input_raw/CE8_S8_L004_R1_001.fastq.gz > input_raw/CE8_S8_L14_R1_001.fastq.gz
cat input_raw/CE8_S8_L001_R2_001.fastq.gz input_raw/CE8_S8_L002_R2_001.fastq.gz input_raw/CE8_S8_L003_R2_001.fastq.gz input_raw/CE8_S8_L004_R2_001.fastq.gz > input_raw/CE8_S8_L14_R2_001.fastq.gz

