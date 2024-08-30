#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



cellranger count --id=H1_hESC_cellranger \
                   --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/H1_hESC \
                   --sample=CE1

