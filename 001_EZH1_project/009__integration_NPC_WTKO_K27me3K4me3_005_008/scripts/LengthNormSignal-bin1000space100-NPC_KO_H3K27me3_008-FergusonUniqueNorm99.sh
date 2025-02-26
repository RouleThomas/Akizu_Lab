#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



computeMatrix scale-regions -S output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99.bw \
    -R ../../Master/meta/GRCh38_bin1000space100.bed \
    --outFileName output/edgeR/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal-bin1000space100-NPC_KO_H3K27me3_008-FergusonUniqueNorm99.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 1000 \
    --regionBodyLength 1000



