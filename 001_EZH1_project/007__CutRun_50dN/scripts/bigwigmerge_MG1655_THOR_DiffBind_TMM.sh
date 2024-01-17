#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
## WT
wiggletools write_bg output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bedGraph \
    median output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1-rep0.bw \
    output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1-rep1.bw
## KO
wiggletools write_bg output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bedGraph \
    median output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2-rep0.bw \
    output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2-rep1.bw
## KOEF1a
wiggletools write_bg output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bedGraph \
    median output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2-rep0.bw \
    output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2-rep1.bw




# Sort the bedgraph 
bedtools sort -i output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bedGraph > \
    output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.sorted.bedGraph
bedtools sort -i output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bedGraph > \
    output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.sorted.bedGraph
bedtools sort -i output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bedGraph > \
    output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.sorted.bedGraph

    
# Convert bedgraph to bigwig
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bw
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bw
bedGraphToBigWig output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bw


