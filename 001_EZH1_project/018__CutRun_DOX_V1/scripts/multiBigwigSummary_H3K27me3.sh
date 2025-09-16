#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b \
    output/bigwig/ESC_WT_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_H3K27me3_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_H3K27me3_R3.unique.dupmark.sorted.bw \
    \
    output/bigwig/ESC_KO_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_KO_H3K27me3_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_KO_H3K27me3_R3.unique.dupmark.sorted.bw \
    \
    output/bigwig/ESC_OEKO_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_OEKO_H3K27me3_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_OEKO_H3K27me3_R3.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R3.unique.dupmark.sorted.bw \
 -o output/bigwig/multiBigwigSummary_H3K27me3.npz

