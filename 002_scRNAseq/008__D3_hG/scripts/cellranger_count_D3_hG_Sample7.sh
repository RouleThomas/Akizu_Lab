#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



                   
cellranger count --id=Sample7 \
  --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
  --fastqs=input/ \
  --sample=SIGAD8


                 
                 