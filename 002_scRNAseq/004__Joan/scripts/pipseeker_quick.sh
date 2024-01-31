#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=10

/scr1/users/roulet/Akizu_Lab/Master/software/pipseeker full --fastq input/. \
    --star-index-path meta/pipseeker-gex-reference-GRCm39-2022.04/ \
    --output-path output_quick/ \
    --skip-version-check \
    --threads 0
