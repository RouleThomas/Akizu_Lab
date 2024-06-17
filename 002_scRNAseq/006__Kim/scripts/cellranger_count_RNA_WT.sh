#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00

cellranger count --id=RNA_WT \
                   --transcriptome=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-gex-mm10-2020-A \
                   --fastqs=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/input_raw/F001/B1 \
                   --chemistry ARC-v1 \
                   --sample=B1_CKDL230036293-1A_22FCTTLT3 # all strings before the `S*_L00*` 
                   




