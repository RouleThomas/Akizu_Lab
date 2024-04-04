#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/NPC_WT_H3K27ac_median.bedGraph \
    median output/bigwig/NPC_WT_H3K27ac_R1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_H3K27ac_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/NPC_WT_H3K27me3_median.bedGraph \
    median output/bigwig/NPC_WT_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_H3K27me3_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/NPC_WT_H3K4me3_median.bedGraph \
    median output/bigwig/NPC_WT_H3K4me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_H3K4me3_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/NPC_WT_H3K9me3_median.bedGraph \
    median output/bigwig/NPC_WT_H3K9me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_H3K9me3_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/NPC_WT_IGG_median.bedGraph \
    median output/bigwig/NPC_WT_IGG_R1.unique.dupmark.sorted.bw \
    output/bigwig/NPC_WT_IGG_R2.unique.dupmark.sorted.bw

# Sort the bedgraph 
bedtools sort -i output/bigwig/NPC_WT_H3K27ac_median.bedGraph > \
    output/bigwig/NPC_WT_H3K27ac_median.sorted.bedGraph

bedtools sort -i output/bigwig/NPC_WT_H3K27me3_median.bedGraph > \
    output/bigwig/NPC_WT_H3K27me3_median.sorted.bedGraph
    
bedtools sort -i output/bigwig/NPC_WT_H3K4me3_median.bedGraph > \
    output/bigwig/NPC_WT_H3K4me3_median.sorted.bedGraph

bedtools sort -i output/bigwig/NPC_WT_H3K9me3_median.bedGraph > \
    output/bigwig/NPC_WT_H3K9me3_median.sorted.bedGraph

bedtools sort -i output/bigwig/NPC_WT_IGG_median.bedGraph > \
    output/bigwig/NPC_WT_IGG_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/NPC_WT_H3K27ac_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/NPC_WT_H3K27ac_median.bw

bedGraphToBigWig output/bigwig/NPC_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/NPC_WT_H3K27me3_median.bw

bedGraphToBigWig output/bigwig/NPC_WT_H3K4me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/NPC_WT_H3K4me3_median.bw

bedGraphToBigWig output/bigwig/NPC_WT_H3K9me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/NPC_WT_H3K9me3_median.bw

bedGraphToBigWig output/bigwig/NPC_WT_IGG_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/NPC_WT_IGG_median.bw


