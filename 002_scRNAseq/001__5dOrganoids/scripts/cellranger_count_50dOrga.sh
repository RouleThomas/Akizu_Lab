#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00



cellranger count --id=count_50dOrga \
                   --transcriptome=meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input_50dOrga/ \
                   --sample=SRR8734991, SRR10914868

                   