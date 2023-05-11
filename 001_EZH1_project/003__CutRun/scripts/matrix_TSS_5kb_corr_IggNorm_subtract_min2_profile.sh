#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_TSS_5kb_corr_IggNorm_subtract_min2.gz \
    -out output/deeptools/matrix_TSS_5kb_corr_IggNorm_subtract_min2_profile.png \
    --perGroup \
    --colors black blue red \
    --plotTitle "" --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


