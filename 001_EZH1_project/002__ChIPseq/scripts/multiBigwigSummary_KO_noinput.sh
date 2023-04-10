#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/2dN_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig/2dN_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig/ESC_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig/ESC_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig/NPC_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig/NPC_KO_H3K27me3_R2.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_KO_noinput.npz


