#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_014R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_EZH2_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_EZH2_014R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_EZH2_014R2_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_EZH2_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_EZH2_010R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_EZH2_014R1_unique_norm99.bw -o output/bigwig_Ferguson/multiBigwigSummary_EZH2_unique_norm99.npz




