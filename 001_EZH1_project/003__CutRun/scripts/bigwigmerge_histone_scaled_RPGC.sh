#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_histone_RPGC/8wN_HET_H3K27me3_median.bedGraph \
    median output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R4.dupmark.sorted.bw
wiggletools write_bg output/bigwig_histone_RPGC/8wN_HET_IGG_median.bedGraph \
    median output/bigwig_histone_RPGC/8wN_HET_IGG_R1.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_HET_IGG_R2.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_HET_IGG_R3.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_HET_IGG_R4.dupmark.sorted.bw
wiggletools write_bg output/bigwig_histone_RPGC/8wN_KO_H3K27me3_median.bedGraph \
    median output/bigwig_histone_RPGC/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_KO_H3K27me3_R4.dupmark.sorted.bw
wiggletools write_bg output/bigwig_histone_RPGC/8wN_KO_IGG_median.bedGraph \
    median output/bigwig_histone_RPGC/8wN_KO_IGG_R1.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_KO_IGG_R2.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_KO_IGG_R3.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_KO_IGG_R4.dupmark.sorted.bw   
wiggletools write_bg output/bigwig_histone_RPGC/8wN_WT_H3K27me3_median.bedGraph \
    median output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R4.dupmark.sorted.bw   
wiggletools write_bg output/bigwig_histone_RPGC/8wN_WT_IGG_median.bedGraph \
    median output/bigwig_histone_RPGC/8wN_WT_IGG_R1.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_WT_IGG_R2.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_WT_IGG_R3.dupmark.sorted.bw \
    output/bigwig_histone_RPGC/8wN_WT_IGG_R4.dupmark.sorted.bw   



# Sort the bedgraph 
bedtools sort -i output/bigwig_histone_RPGC/8wN_HET_H3K27me3_median.bedGraph > \
    output/bigwig_histone_RPGC/8wN_HET_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_histone_RPGC/8wN_HET_IGG_median.bedGraph > \
    output/bigwig_histone_RPGC/8wN_HET_IGG_median.sorted.bedGraph

bedtools sort -i output/bigwig_histone_RPGC/8wN_KO_H3K27me3_median.bedGraph > \
    output/bigwig_histone_RPGC/8wN_KO_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_histone_RPGC/8wN_KO_IGG_median.bedGraph > \
    output/bigwig_histone_RPGC/8wN_KO_IGG_median.sorted.bedGraph

bedtools sort -i output/bigwig_histone_RPGC/8wN_WT_H3K27me3_median.bedGraph > \
    output/bigwig_histone_RPGC/8wN_WT_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_histone_RPGC/8wN_WT_IGG_median.bedGraph > \
    output/bigwig_histone_RPGC/8wN_WT_IGG_median.sorted.bedGraph


    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_histone_RPGC/8wN_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_RPGC/8wN_HET_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_histone_RPGC/8wN_HET_IGG_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_RPGC/8wN_HET_IGG_median.bw

bedGraphToBigWig output/bigwig_histone_RPGC/8wN_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_RPGC/8wN_KO_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_histone_RPGC/8wN_KO_IGG_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_RPGC/8wN_KO_IGG_median.bw

bedGraphToBigWig output/bigwig_histone_RPGC/8wN_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_RPGC/8wN_WT_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_histone_RPGC/8wN_WT_IGG_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_RPGC/8wN_WT_IGG_median.bw