#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/2dN_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig/2dN_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig/ESC_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig/ESC_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig/NPC_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig/NPC_HET_H3K27me3_R2.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_HET_noinput.npz


