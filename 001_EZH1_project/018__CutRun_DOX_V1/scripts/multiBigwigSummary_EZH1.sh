#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b \
    output/bigwig/ESC_WT_EZH1_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_EZH1_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_EZH1_R3.unique.dupmark.sorted.bw \
    \
    output/bigwig/ESC_KO_EZH1_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_KO_EZH1_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_KO_EZH1_R3.unique.dupmark.sorted.bw \
    \
    output/bigwig/ESC_OEKO_EZH1_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_OEKO_EZH1_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_OEKO_EZH1_R3.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R1.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R2.unique.dupmark.sorted.bw \
    output/bigwig/ESC_WT_IGG_R3.unique.dupmark.sorted.bw \
 -o output/bigwig/multiBigwigSummary_EZH1.npz

