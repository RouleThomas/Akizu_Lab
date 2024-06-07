#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00

cellranger count --id=Kcnc1_p180_CX_Rep1 \
                   --transcriptome=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-gex-mm10-2020-A \
                   --fastqs=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/005__Goldberg/input_raw/snRNAseq_Kcnc1_R320H/snRNAseq_Kcnc1_reorganized/Kcnc1_P14/Kcnc1_P180_CTX_Sample1 \
                   --sample=KCNC1-SN-P180-M1-CTX-EthanGoldberg-02022023  # all strings before the `S*_L00*`




