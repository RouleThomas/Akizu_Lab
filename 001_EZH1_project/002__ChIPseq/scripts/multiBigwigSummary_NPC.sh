#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b  output/bigwig/NPC_HET_input_R1.dupmark.sorted.bw output/bigwig/NPC_HET_input_R2.dupmark.sorted.bw output/bigwig/NPC_KO_input_R1.dupmark.sorted.bw output/bigwig/NPC_KO_input_R2.dupmark.sorted.bw output/bigwig/NPC_WT_input_R1.dupmark.sorted.bw output/bigwig/NPC_WT_input_R2.dupmark.sorted.bw output/bigwig/NPC_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig/NPC_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig/NPC_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig/NPC_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_R1.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_R2.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_NPC.npz


