#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=100:00:00



computeMatrix scale-regions -S output/bigwig_Ferguson/NPC_KO_H3K27me3_008_unique_norm99_initialBigwig.bw \
    -R output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K27me3_pool_peaks.sorted.merge.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99initialBigwig.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99initialBigwig.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_H3K27me3_008-FergusonUniqueNorm99initialBigwig.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 100 \
    --regionBodyLength 100



