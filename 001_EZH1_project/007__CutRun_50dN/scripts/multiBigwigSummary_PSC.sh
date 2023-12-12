#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig/PSC_WT_EZH1cs_01FA.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_EZH1cs_1FA.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_H3K27me1_01FA.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_H3K27me1_1FA.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_PSC.npz








