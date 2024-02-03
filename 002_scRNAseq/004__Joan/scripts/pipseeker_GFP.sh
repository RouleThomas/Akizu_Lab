#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=10

/scr1/users/roulet/Akizu_Lab/Master/software/pipseeker full --fastq input/. \
    --star-index-path meta/STAR2.7.10b_mm10-2020-A_tdT_EGFP__uniquedl/ \
    --output-path output_EGFP/ \
    --skip-version-check \
    --threads 0
    

