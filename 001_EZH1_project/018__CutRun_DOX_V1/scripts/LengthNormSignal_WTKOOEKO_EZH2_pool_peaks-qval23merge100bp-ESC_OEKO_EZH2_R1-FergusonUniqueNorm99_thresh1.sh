#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



computeMatrix scale-regions -S output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_unique_norm99_initialBigwig_thresh1.bw \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_EZH2_pool_peaks.sorted.merge100bp.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99_thresh1.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99_thresh1.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKOOEKO_EZH2_pool_peaks-qval23merge100bp-ESC_OEKO_EZH2_R1-FergusonUniqueNorm99_thresh1.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 100 \
    --regionBodyLength 100



