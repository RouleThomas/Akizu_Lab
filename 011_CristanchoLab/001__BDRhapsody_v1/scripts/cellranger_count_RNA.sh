#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00

cellranger count --id=RNA_cellranger \
                   --transcriptome=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-gex-mm10-2020-A \
                   --fastqs=/scr1/users/roulet/Akizu_Lab/011_CristanchoLab/001__BDRhapsody_v1/input \
                   --chemistry ARC-v1 \
                   --sample=AC_WTA_SMK_index_library # all strings before the `S*_L00*` 
                   




