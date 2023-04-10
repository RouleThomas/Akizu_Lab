#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/2dN_HET_input_R1.dupmark.sorted.bw output/bigwig/2dN_HET_input_R2.dupmark.sorted.bw output/bigwig/2dN_KO_input_R1.dupmark.sorted.bw output/bigwig/2dN_KO_input_R2.dupmark.sorted.bw output/bigwig/2dN_WT_input_R1.dupmark.sorted.bw output/bigwig/2dN_WT_input_R2.dupmark.sorted.bw output/bigwig/2dN_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig/2dN_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig/2dN_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig/2dN_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig/2dN_WT_H3K27me3_R1.dupmark.sorted.bw output/bigwig/2dN_WT_H3K27me3_R2.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_2dN.npz


