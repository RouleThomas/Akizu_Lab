#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



cellranger count --id=72hgastruloidhumanXMU \
                   --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/72hgastruloidhumanXMU \
                   --sample=72h_XMU

                   