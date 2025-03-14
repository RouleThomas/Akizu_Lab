#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



# Convert to bedGraph and smooth
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_EZH2_006R_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_EZH2_006R_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_EZH2_006R_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_EZH2_010R_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_EZH2_010R_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_EZH2_010R_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_EZH2_014R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_EZH2_014R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_EZH2_014R1_unique_norm99_smooth50bp.bedGraph



multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_EZH2_013R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KO_EZH2_013R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_EZH2_013R1_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_EZH2_014R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KO_EZH2_014R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_EZH2_014R1_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_EZH2_014R2_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KO_EZH2_014R2_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_EZH2_014R2_unique_norm99_smooth50bp.bedGraph


multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_006R_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_006R_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_006R_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_013R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_013R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_013R1_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_014R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_014R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_014R1_unique_norm99_smooth50bp.bedGraph


