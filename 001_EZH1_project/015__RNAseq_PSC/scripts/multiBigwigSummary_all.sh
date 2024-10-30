#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b  output/bigwig/PSC_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_WT_R3_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KO_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KO_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KO_R3_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KOEF1aEZH1_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KOEF1aEZH1_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KOEF1aEZH1_R3_Aligned.sortedByCoord.out.bw -o output/bigwig/multiBigwigSummary_all.npz





