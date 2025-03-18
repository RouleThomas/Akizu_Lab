#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_initialBigwig_median.bedGraph \
    median output/bigwig_Ferguson/NPC_WT_H3K4me3_008_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/NPC_WT_H3K4me3_005_unique_norm99_initialBigwig.bw 
wiggletools write_bg output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_initialBigwig_median.bedGraph \
    median output/bigwig_Ferguson/NPC_KO_H3K4me3_008_unique_norm99_initialBigwig.bw \
    output/bigwig_Ferguson/NPC_KO_H3K4me3_005_unique_norm99_initialBigwig.bw 


# Sort the bedgraph 
bedtools sort -i output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_initialBigwig_median.bedGraph > output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_initialBigwig_median.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_initialBigwig_median.bedGraph > output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_initialBigwig_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_initialBigwig_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_initialBigwig_median.bw

bedGraphToBigWig output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_initialBigwig_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_initialBigwig_median.bw




