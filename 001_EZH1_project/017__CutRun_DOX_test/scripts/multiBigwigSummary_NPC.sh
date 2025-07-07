#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig/NPC_KO_EZH1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_KO_EZH2.unique.dupmark.sorted.bw \
    output/bigwig/NPC_KO_H3K27me3.unique.dupmark.sorted.bw \
    output/bigwig/NPC_KO_IGG.unique.dupmark.sorted.bw \
    output/bigwig/NPC_OEKO_EZH1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_OEKO_EZH2.unique.dupmark.sorted.bw \
    output/bigwig/NPC_OEKO_H3K27me3.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_EZH1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_EZH2.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_H3K27me3.unique.dupmark.sorted.bw \
 -o output/bigwig/multiBigwigSummary_NPC.npz



