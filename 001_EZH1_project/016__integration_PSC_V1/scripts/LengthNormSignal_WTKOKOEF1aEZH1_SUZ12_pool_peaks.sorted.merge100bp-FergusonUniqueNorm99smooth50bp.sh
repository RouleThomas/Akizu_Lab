#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


computeMatrix scale-regions -S output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth50bp.bw \
    -R output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp.bed \
    --outFileName output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp-FergusonUniqueNorm99smooth50bp.gz \
    --outFileNameMatrix output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp-FergusonUniqueNorm99smooth50bp.txt \
    --outFileSortedRegions output/edgeR/LengthNormSignal_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge100bp-FergusonUniqueNorm99smooth50bp.bed













