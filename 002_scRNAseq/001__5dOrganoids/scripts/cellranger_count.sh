#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00



cellranger count --id=count \
                   --transcriptome=meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/ \
                   --sample=SRR8734990

                   