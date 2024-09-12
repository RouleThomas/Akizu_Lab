#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/PSC_WT_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_H3K27me3_R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_H3K27me3_R2.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_WTvsKO_H3K27me3.npz





