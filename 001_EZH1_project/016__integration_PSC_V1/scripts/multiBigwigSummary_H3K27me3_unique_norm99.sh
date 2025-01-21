#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99.bw -o output/bigwig_Ferguson/multiBigwigSummary_H3K27me3_unique_norm99.npz





