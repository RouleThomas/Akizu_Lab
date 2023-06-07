#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00



cellranger count --id=count_SRR10914868 \
                   --transcriptome=meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input_50dOrga/ \
                   --sample=SRR10914868
                   