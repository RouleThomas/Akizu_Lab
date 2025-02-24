#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00



# Convert back bedGraph to bigwig
bedtools sort -i output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_median_smooth50bp.bedGraph > output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_median_smooth50bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_median_smooth50bp.bedGraph > output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_median_smooth50bp.sorted.bedGraph

bedtools sort -i output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median_smooth50bp.bedGraph > output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median_smooth50bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median_smooth50bp.bedGraph > output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median_smooth50bp.sorted.bedGraph


bedGraphToBigWig \
    output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_median_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/NPC_WT_H3K4me3_unique_norm99_median_smooth50bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_median_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/NPC_KO_H3K4me3_unique_norm99_median_smooth50bp.bw

bedGraphToBigWig \
    output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median_smooth50bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median_smooth50bp.bw

