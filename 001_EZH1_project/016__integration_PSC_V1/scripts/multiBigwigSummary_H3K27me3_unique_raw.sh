#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bw \
    output/bigwig/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bw \
    output/bigwig/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_H3K27me3_unique_raw.npz





