#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8



multiBigwigSummary bins -b output/bigwig/50dN_WT_EZH1.unique.dupmark.sorted.bw \
output/bigwig/50dN_WT_EZH2.unique.dupmark.sorted.bw \
output/bigwig/50dN_WT_H3K27ac.unique.dupmark.sorted.bw \
output/bigwig/50dN_WT_H3K27me1AM.unique.dupmark.sorted.bw \
output/bigwig/50dN_WT_H3K27me1OR.unique.dupmark.sorted.bw \
output/bigwig/50dN_WT_H3K27me3.unique.dupmark.sorted.bw \
output/bigwig/50dN_WT_IGG.unique.dupmark.sorted.bw \
output/bigwig/50dN_WT_SUZ12.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_all.npz -p max








