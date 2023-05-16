#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_TSS_5kb_WT_DiffBind_TMM.gz \
    -out output/deeptools/matrix_TSS_5kb_WT_DiffBind_TMM_profile.png \
    --perGroup \
    --plotTitle "" --samplesLabel "Rep1" "Rep2" "Rep3" "Rep4" \
    --refPointLabel "TSS" \
    -T "WT read density" \
    -z ""


