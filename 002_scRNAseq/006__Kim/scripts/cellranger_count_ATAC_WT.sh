#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00

cellranger-atac count --id=ATAC_WT \
                   --reference=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                   --fastqs=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/input_raw/F002/A1 \
                   --chemistry ARC-v1 \
                   --sample=A1SubLib_CKDL230036285-1A_22FF27LT3  # all strings before the `S*_L00*`




