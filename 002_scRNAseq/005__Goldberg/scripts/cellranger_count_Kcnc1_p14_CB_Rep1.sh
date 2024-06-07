#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00

cellranger count --id=Kcnc1_p14_CB_Rep1 \
                   --transcriptome=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-gex-mm10-2020-A \
                   --fastqs=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/005__Goldberg/input_raw/snRNAseq_Kcnc1_R320H/8_31_23_core_reanalyzed_corrupted_files/fastq/Kcnc1_P14_CERE_Sample1 \
                   --sample=KCNC1-SN-P14-M1-1-CER  # all strings before the `S*_L00*`




