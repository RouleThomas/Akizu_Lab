#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_005R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_SUZ12_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_SUZ12_014R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_SUZ12_014R2_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_SUZ12_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_SUZ12_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_SUZ12_014R1_unique_norm99.bw -o output/bigwig_Ferguson/multiBigwigSummary_SUZ12_unique_norm99.npz



