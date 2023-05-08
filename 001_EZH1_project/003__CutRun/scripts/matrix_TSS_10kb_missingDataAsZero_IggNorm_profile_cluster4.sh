#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00





plotProfile -m output/deeptools/matrix_TSS_10kb_missingDataAsZero_IggNorm.gz \
    -out output/deeptools/matrix_TSS_10kb_missingDataAsZero_IggNorm_cluster4.png \
    --kmeans 4 \
    --numPlotsPerRow 4 \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



