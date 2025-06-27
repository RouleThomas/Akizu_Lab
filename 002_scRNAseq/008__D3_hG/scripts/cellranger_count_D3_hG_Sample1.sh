#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



                   
cellranger count --id=Sample1 \
  --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
  --fastqs=input/ \
  --sample=SIGAA8


