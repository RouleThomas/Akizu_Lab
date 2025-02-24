#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00





# Convert to bedGraph and smooth

multiBigwigSummary bins -b output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_median.bw --binSize 50 -o output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_median_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_median_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_median.bw --binSize 50 -o output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_median_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_median_smooth50bp.bedGraph

multiBigwigSummary bins -b output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median.bw --binSize 50 -o output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median.bw --binSize 50 -o output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median_smooth50bp.bedGraph










