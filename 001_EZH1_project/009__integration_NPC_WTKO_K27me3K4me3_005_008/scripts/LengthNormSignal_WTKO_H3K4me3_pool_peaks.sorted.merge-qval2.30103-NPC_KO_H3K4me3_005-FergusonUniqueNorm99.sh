#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=100:00:00



computeMatrix scale-regions -S output/bigwig_Ferguson/NPC_KO_H3K4me3_005_unique_norm99.bw \
    -R output/macs2/broad/broad_blacklist_qval2.30103/NPC_WTKO_H3K4me3_pool_peaks.sorted.merge.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKO_H3K4me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K4me3_005-FergusonUniqueNorm99.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKO_H3K4me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K4me3_005-FergusonUniqueNorm99.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKO_H3K4me3_pool_peaks.sorted.merge-qval2.30103-NPC_KO_H3K4me3_005-FergusonUniqueNorm99.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 100 \
    --regionBodyLength 100



