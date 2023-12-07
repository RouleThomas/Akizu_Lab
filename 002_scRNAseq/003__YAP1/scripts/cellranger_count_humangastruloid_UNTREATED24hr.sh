#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



cellranger count --id=24hgastruloidhumanUN \
                   --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/24hgastruloidhumanUN \
                   --sample=24hgastruloidhumanUN

                   