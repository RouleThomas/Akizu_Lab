#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/ESC_HET_input_R1.dupmark.sorted.bw output/bigwig/ESC_HET_input_R2.dupmark.sorted.bw output/bigwig/ESC_KO_input_R1.dupmark.sorted.bw output/bigwig/ESC_KO_input_R2.dupmark.sorted.bw output/bigwig/ESC_WT_input_R1.dupmark.sorted.bw output/bigwig/ESC_WT_input_R2.dupmark.sorted.bw output/bigwig/ESC_WT_input_R3.dupmark.sorted.bw output/bigwig/ESC_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig/ESC_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig/ESC_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig/ESC_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig/ESC_WT_H3K27me3_R1.dupmark.sorted.bw output/bigwig/ESC_WT_H3K27me3_R2.dupmark.sorted.bw output/bigwig/ESC_WT_H3K27me3_R3.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_ESC.npz


