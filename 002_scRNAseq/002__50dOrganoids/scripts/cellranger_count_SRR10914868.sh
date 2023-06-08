#!/bin/bash
#SBATCH --mem=750G
#SBATCH --time=100:00:00



/scr1/users/roulet/Akizu_Lab/Master/software/cellranger-6.0.2/cellranger count --id=count_SRR10914868 \
                   --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/ \
                   --sample=SRR10914868
                   