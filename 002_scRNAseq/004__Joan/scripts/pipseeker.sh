#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


/scr1/users/roulet/Akizu_Lab/Master/software/pipseeker full --fastq input/. \
    --star-index-path meta/pipseeker-gex-reference-GRCm39-2022.04/ \
    --output-path output/ \
    --skip-version-check 
