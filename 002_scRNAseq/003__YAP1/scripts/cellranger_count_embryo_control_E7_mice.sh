#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=100:00:00



cellranger count --id=E7mousecontrol \
                   --transcriptome=../meta/refdata-gex-mm10-2020-A \
                   --fastqs=input/E7mousecontrol \
                   --sample=E7mousecontrol

                   