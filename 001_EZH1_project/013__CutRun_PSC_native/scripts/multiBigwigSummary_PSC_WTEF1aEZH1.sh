#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/PSC_WTEF1aEZH1_EZH1_R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_EZH1_R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_EZH2_R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_EZH2_R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_H3K27me3_R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_IGG_R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_IGG_R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_SUZ12_R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WTEF1aEZH1_SUZ12_R2.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_PSC_WTEF1aEZH1.npz








