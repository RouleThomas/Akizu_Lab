#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig/ESC_KO_EZH1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_KO_EZH2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_KO_H3K27me3.unique.dupmark.sorted.bw \
    output/bigwig/ESC_KO_H3K27me3epi.unique.dupmark.sorted.bw \
    output/bigwig/ESC_OEKO_EZH1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_OEKO_EZH2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_OEKO_H3K27me3.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_EZH1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_EZH2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_H3K27me3.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_H3K27me3epi.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG.unique.dupmark.sorted.bw \
 -o output/bigwig/multiBigwigSummary_ESC.npz



