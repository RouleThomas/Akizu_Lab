#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00





multiBigwigSummary bins -b output/bigwig/PSC_KOEF1aEZH1_SUZ12.unique.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_EZH2.unique.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_HA.unique.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_IGG.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_SUZ12.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_EZH2.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_HA.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_H3K27me3.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_IGG.unique.dupmark.sorted.bw \
output/bigwig/PSC_KO_SUZ12.unique.dupmark.sorted.bw \
output/bigwig/PSC_KO_EZH2.unique.dupmark.sorted.bw \
output/bigwig/PSC_KO_HA.unique.dupmark.sorted.bw \
output/bigwig/PSC_KO_EZH1cs.unique.dupmark.sorted.bw \
output/bigwig/PSC_KO_H3K27me3.unique.dupmark.sorted.bw \
output/bigwig/PSC_KO_IGG.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_all.npz








