#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=8



cellranger-arc count --id=multiome_Bap1KO \
                       --reference=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/meta/refdata-cellranger-arc-mm10-2020-A-2.0.0 \
                       --libraries=/scr1/users/roulet/Akizu_Lab/002_scRNAseq/006__Kim/libraries_Bap1KO.csv \
                       --localcores=8 \
                       --localmem=500


                       