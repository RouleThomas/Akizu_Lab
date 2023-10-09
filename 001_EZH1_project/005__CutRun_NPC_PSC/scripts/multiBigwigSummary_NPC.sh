#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig/NPC_WT_EZH1cs.dupmark.sorted.bw \
output/bigwig/NPC_KO_EZH1cs.dupmark.sorted.bw \
output/bigwig/NPC_WT_EZH1pt.dupmark.sorted.bw \
output/bigwig/NPC_KO_EZH1pt.dupmark.sorted.bw \
output/bigwig/NPC_WT_H3K27me1.dupmark.sorted.bw \
output/bigwig/NPC_KO_H3K27me1.dupmark.sorted.bw \
output/bigwig/NPC_WT_H3K27me3.dupmark.sorted.bw \
output/bigwig/NPC_KO_H3K27me3.dupmark.sorted.bw \
output/bigwig/NPC_WT_H3K4me3.dupmark.sorted.bw \
output/bigwig/NPC_KO_H3K4me3.dupmark.sorted.bw \
output/bigwig/NPC_WT_EZH2.dupmark.sorted.bw \
output/bigwig/NPC_KO_EZH2.dupmark.sorted.bw \
output/bigwig/NPC_WT_SUZ12.dupmark.sorted.bw \
output/bigwig/NPC_KO_SUZ12.dupmark.sorted.bw \
output/bigwig/NPC_WT_IGG.dupmark.sorted.bw \
output/bigwig/NPC_KO_IGG.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_NPC.npz








