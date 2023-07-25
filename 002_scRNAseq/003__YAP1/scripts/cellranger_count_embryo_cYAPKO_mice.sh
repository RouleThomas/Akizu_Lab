#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=100:00:00



cellranger count --id=embryo_cYAPKO_e775_mice \
                   --transcriptome=../meta/refdata-gex-mm10-2020-A \
                   --fastqs=input/embryo_cYAPKO_e775_L1-ds.1b49bd37e60f41349c1ed68a0f5c130e \
                   --sample=embryo_cYAPKO_e775

                   