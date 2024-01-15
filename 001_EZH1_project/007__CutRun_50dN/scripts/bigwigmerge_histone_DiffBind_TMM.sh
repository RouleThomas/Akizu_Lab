#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_histone_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_median.bedGraph \
    median output/bigwig_histone_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig_histone_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig_histone_DiffBind_TMM/50dN_KO_H3K27me3_median.bedGraph \
    median output/bigwig_histone_DiffBind_TMM/50dN_KO_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig_histone_DiffBind_TMM/50dN_KO_H3K27me3_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_median.bedGraph \
    median output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bw




# Sort the bedgraph 
bedtools sort -i output/bigwig_histone_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_median.bedGraph > \
    output/bigwig_histone_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_histone_DiffBind_TMM/50dN_KO_H3K27me3_median.bedGraph > \
    output/bigwig_histone_DiffBind_TMM/50dN_KO_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_median.bedGraph > \
    output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_median.sorted.bedGraph

    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_histone_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_histone_DiffBind_TMM/50dN_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_DiffBind_TMM/50dN_KO_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_median.bw


