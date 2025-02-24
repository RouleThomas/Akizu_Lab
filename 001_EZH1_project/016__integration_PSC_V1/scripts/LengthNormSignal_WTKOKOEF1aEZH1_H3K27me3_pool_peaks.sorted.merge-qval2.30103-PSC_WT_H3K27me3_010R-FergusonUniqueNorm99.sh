#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



computeMatrix scale-regions -S output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bw \
    -R output/macs2/broad/broad_blacklist_qval2.30103/PSC_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_H3K27me3_pool_peaks.sorted.merge-qval2.30103-PSC_WT_H3K27me3_010R-FergusonUniqueNorm99.bed \
    --missingDataAsZero \
    --averageTypeBins sum \
    --binSize 100 \
    --regionBodyLength 100



