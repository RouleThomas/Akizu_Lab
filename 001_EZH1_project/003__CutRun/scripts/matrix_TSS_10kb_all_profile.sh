#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_TSS_10kb_all.gz \
    -out output/deeptools/matrix_TSS_10kb_all_profile.png \
    --perGroup \
    --colors black blue red darkred \
    --plotTitle "" --samplesLabel "WT" "HET" "KO" "patient" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



