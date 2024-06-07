#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00

cellranger count --id=WT_p35_CX_Rep1 \
                   --transcriptome=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-gex-mm10-2020-A \
                   --fastqs=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/005__Goldberg/input_raw/snRNAseq_Kcnc1_R320H/8_31_23_core_reanalyzed_corrupted_files/fastq/WT_P35_CTX_Sample1 \
                   --sample=KCNC1-SN-P35-W1-CTX-EthanGoldberg-01232023  # all strings before the `S*_L00*`




