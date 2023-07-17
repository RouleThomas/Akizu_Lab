#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=100:00:00



cellranger count --id=embryo_control_e775 \
                   --transcriptome=../meta/refdata-gex-GRCh38-2020-A \
                   --fastqs=input/embryo_control_e775_L1-ds.7627bbbce9dd4313826228c04bae76eb \
                   --sample=embryo_control_e775

                   