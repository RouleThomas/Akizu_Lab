#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



multiBigwigSummary bins -b output/bigwig_downsample/2dN_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/2dN_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/2dN_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/2dN_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/2dN_WT_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/2dN_WT_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/ESC_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/ESC_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/ESC_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/ESC_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/ESC_WT_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/ESC_WT_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/ESC_WT_H3K27me3_R3.dupmark.sorted.bw output/bigwig_downsample/NPC_HET_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/NPC_HET_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/NPC_KO_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/NPC_KO_H3K27me3_R2.dupmark.sorted.bw output/bigwig_downsample/NPC_WT_H3K27me3_R1.dupmark.sorted.bw output/bigwig_downsample/NPC_WT_H3K27me3_R2.dupmark.sorted.bw -o output/bigwig_downsample/multibigwigSummary_all_H3K27me3_downsample.npz

