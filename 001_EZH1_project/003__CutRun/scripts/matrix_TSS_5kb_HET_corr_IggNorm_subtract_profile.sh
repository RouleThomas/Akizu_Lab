#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_TSS_5kb_HET_corr_IggNorm_subtract.gz \
    -out output/deeptools/matrix_TSS_5kb_HET_corr_IggNorm_subtract_profile.png \
    --perGroup \
    --plotTitle "" --samplesLabel "Rep1" "Rep2" "Rep3" "Rep4" \
    --refPointLabel "TSS" \
    -T "HET read density" \
    -z ""


