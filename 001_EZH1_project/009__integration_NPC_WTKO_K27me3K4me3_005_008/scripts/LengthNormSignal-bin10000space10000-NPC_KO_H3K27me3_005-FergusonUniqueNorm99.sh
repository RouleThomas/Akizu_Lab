#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



computeMatrix scale-regions -S output/bigwig_Ferguson/NPC_KO_H3K27me3_005_unique_norm99.bw \
    -R ../../Master/meta/GRCh38_bin10000space10000.bed \
    --outFileName output/edgeR/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal-bin10000space10000-NPC_KO_H3K27me3_005-FergusonUniqueNorm99.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 10000 \
    --regionBodyLength 10000



