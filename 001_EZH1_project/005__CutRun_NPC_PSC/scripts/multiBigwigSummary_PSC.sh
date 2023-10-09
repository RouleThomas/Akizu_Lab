#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


multiBigwigSummary bins -b output/bigwig/PSC_KOEF1aEZH1_EZH1cs.dupmark.sorted.bw \
output/bigwig/PSC_KOsynEZH1_EZH1cs.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_EZH1pt.dupmark.sorted.bw \
output/bigwig/PSC_KOsynEZH1_EZH1pt.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_H3K27me3.dupmark.sorted.bw \
output/bigwig/PSC_KOsynEZH1_H3K27me3.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_HA.dupmark.sorted.bw \
output/bigwig/PSC_KOsynEZH1_HA.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_SUZ12.dupmark.sorted.bw \
output/bigwig/PSC_KOsynEZH1_SUZ12.dupmark.sorted.bw \
output/bigwig/PSC_KOEF1aEZH1_IGG.dupmark.sorted.bw \
output/bigwig/PSC_KOsynEZH1_IGG.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_PSC.npz








