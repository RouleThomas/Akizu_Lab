#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00




multiBigwigSummary bins -b output/bigwig/50dN_KOEF1aEZH1_EZH1cs_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_EZH1cs_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_EZH2_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_EZH2_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_H3K27me3_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_H3K27me3_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_IGG_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_SUZ12_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KOEF1aEZH1_SUZ12_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_EZH1cs_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_EZH1cs_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_EZH2_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_EZH2_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_H3K27me3_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_H3K27me3_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_IGG_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_IGG_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_SUZ12_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_KO_SUZ12_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_EZH1cs_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_EZH1cs_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_EZH2_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_EZH2_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_H3K27me3_R3.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_IGG_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_IGG_R2.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_SUZ12_R1.unique.dupmark.sorted.bw \
output/bigwig/50dN_WTQ731E_SUZ12_R2.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_EZH1cs_01FA.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_EZH1cs_1FA.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_H3K27me1_01FA.unique.dupmark.sorted.bw \
output/bigwig/PSC_WT_H3K27me1_1FA.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_all.npz








