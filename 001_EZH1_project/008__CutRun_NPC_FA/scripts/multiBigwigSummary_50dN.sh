#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8



multiBigwigSummary bins -b output/bigwig/50dNFA_KOEF1aEZH1_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/50dNnative_KOEF1aEZH1_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/50dNFA_KOEF1aEZH1_EZH2.unique.dupmark.sorted.bw \
output/bigwig/50dNnative_KOEF1aEZH1_EZH2.unique.dupmark.sorted.bw \
output/bigwig/50dNFA_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_50dN.npz -p max








