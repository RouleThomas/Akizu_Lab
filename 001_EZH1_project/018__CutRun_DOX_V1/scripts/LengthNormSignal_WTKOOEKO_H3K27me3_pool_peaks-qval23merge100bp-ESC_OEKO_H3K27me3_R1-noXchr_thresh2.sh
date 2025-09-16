#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00





computeMatrix scale-regions -S output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bw \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_noXchr_pool_peaks.sorted.merge100bp.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-noXchr_thresh2.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-noXchr_thresh2.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKOOEKO_H3K27me3_pool_peaks-qval23merge100bp-ESC_OEKO_H3K27me3_R1-noXchr_thresh2.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 100 \
    --regionBodyLength 100



