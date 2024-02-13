#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


multiBigwigSummary bins -b output/bigwig/NPC_KO_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/NPC_KO_EZH2.unique.dupmark.sorted.bw \
output/bigwig/NPC_KO_H3K27ac.unique.dupmark.sorted.bw \
output/bigwig/NPC_KO_IGG.unique.dupmark.sorted.bw \
output/bigwig/NPC_KOEF1aEZH1_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/NPC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw \
output/bigwig/NPC_KOEF1aEZH1_H3K4me3.unique.dupmark.sorted.bw \
output/bigwig/NPC_KOEF1aEZH1_IGG.unique.dupmark.sorted.bw \
output/bigwig/NPC_KOEF1aEZH1_SUZ12.unique.dupmark.sorted.bw \
output/bigwig/NPC_WT_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/NPC_WT_EZH2.unique.dupmark.sorted.bw \
output/bigwig/NPC_WT_H3K27ac.unique.dupmark.sorted.bw \
output/bigwig/NPC_WT_H3K4me3.unique.dupmark.sorted.bw \
output/bigwig/NPC_WT_SUZ12.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_NPC.npz -p max




