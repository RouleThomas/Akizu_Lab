#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8



multiBigwigSummary bins -b output/bigwig/hESC_WT_QSER1FLAG_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_QSER1FLAG_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_DNMT3A_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_DNMT3A_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_DNMT3B_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_DNMT3B_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_H3K27me3_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_H3K27me3_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_H3K4me3_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_H3K4me3_R2.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_EZH2_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_input_R1.unique.dupmark.sorted.bw \
output/bigwig/hESC_WT_input_R2.unique.dupmark.sorted.bw -o output/bigwig/multiBigwigSummary_all.npz -p max



