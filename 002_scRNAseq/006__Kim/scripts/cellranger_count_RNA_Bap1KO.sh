#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00

cellranger count --id=RNA_Bap1KO \
                   --transcriptome=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-gex-mm10-2020-A \
                   --fastqs=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/input_raw/F001/B2 \
                   --chemistry ARC-v1 \
                   --sample=B2_CKDL230036294-1A_22FCTTLT3  # all strings before the `S*_L00*`




