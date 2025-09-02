#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



computeMatrix scale-regions -S output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_unique_norm99_initialBigwig.bw \
    -R output/macs2/broad/broad_blacklist_qval3/ESC_WTKOOEKO_H3K27me3_pool_peaks.sorted.merge100bp.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval3merge100bp-ESC_WT_H3K27me3_R1-FergusonUniqueNorm99.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 100 \
    --regionBodyLength 100



