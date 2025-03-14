#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00





# Convert to bedGraph and smooth
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth250bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth250bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth250bp.bedGraph



multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth250bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth250bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth250bp.bedGraph


multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth250bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth250bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median.bw --binSize 250 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth250bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth250bp.bedGraph











