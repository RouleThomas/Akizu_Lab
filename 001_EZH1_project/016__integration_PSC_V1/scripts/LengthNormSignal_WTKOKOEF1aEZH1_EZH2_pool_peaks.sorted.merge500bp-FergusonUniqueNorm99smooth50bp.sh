#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


computeMatrix scale-regions -S output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth50bp.bw \
    -R output/macs2/broad/PSC_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp-FergusonUniqueNorm99smooth50bp.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp-FergusonUniqueNorm99smooth50bp.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_EZH2_pool_peaks.sorted.merge500bp-FergusonUniqueNorm99smooth50bp.bed













