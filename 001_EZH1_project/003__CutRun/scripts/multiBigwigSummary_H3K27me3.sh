#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
output/bigwig/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
output/bigwig/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
output/bigwig/8wN_WT_H3K27me3_R4.dupmark.sorted.bw \
output/bigwig/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
output/bigwig/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
output/bigwig/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
output/bigwig/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
output/bigwig/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
output/bigwig/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
output/bigwig/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
output/bigwig/8wN_HET_H3K27me3_R4.dupmark.sorted.bw \
output/bigwig/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bw \
output/bigwig/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_H3K27me3.npz








