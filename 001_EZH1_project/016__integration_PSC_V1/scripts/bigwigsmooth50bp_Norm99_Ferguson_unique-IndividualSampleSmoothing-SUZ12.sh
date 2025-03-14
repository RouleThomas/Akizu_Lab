#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



# Convert to bedGraph and smooth
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_SUZ12_006R_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_SUZ12_006R_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_SUZ12_006R_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_SUZ12_013R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_SUZ12_013R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_SUZ12_013R1_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_SUZ12_014R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_SUZ12_014R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_SUZ12_014R1_unique_norm99_smooth50bp.bedGraph



multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_SUZ12_013R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KO_SUZ12_013R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_SUZ12_013R1_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_SUZ12_014R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KO_SUZ12_014R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_SUZ12_014R1_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KO_SUZ12_014R2_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KO_SUZ12_014R2_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KO_SUZ12_014R2_unique_norm99_smooth50bp.bedGraph



multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_005R_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_005R_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_005R_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_006R_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_006R_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_006R_unique_norm99_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_013R1_unique_norm99.bw --binSize 50 -o output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_013R1_unique_norm99_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_013R1_unique_norm99_smooth50bp.bedGraph


