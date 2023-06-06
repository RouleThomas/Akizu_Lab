#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_DiffBind_LIB/NPC_HET_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/NPC_HET_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/NPC_HET_H3K27me3_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/NPC_KO_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/NPC_KO_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/NPC_KO_H3K27me3_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/NPC_WT_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/NPC_WT_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/NPC_WT_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_DiffBind_LIB/ESC_HET_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/ESC_HET_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/ESC_HET_H3K27me3_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/ESC_KO_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/ESC_KO_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/ESC_KO_H3K27me3_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/ESC_WT_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/ESC_WT_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/ESC_WT_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_DiffBind_LIB/2dN_HET_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/2dN_HET_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/2dN_HET_H3K27me3_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/2dN_KO_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/2dN_KO_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/2dN_KO_H3K27me3_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/2dN_WT_H3K27me3_median.bedGraph \
    median output/bigwig_DiffBind_LIB/2dN_WT_H3K27me3_R1.bw \
    output/bigwig_DiffBind_LIB/2dN_WT_H3K27me3_R2.bw



wiggletools write_bg output/bigwig_DiffBind_LIB/NPC_HET_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/NPC_HET_input_R1.bw \
    output/bigwig_DiffBind_LIB/NPC_HET_input_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/NPC_KO_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/NPC_KO_input_R1.bw \
    output/bigwig_DiffBind_LIB/NPC_KO_input_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/NPC_WT_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/NPC_WT_input_R1.bw \
    output/bigwig_DiffBind_LIB/NPC_WT_input_R2.bw

wiggletools write_bg output/bigwig_DiffBind_LIB/ESC_HET_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/ESC_HET_input_R1.bw \
    output/bigwig_DiffBind_LIB/ESC_HET_input_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/ESC_KO_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/ESC_KO_input_R1.bw \
    output/bigwig_DiffBind_LIB/ESC_KO_input_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/ESC_WT_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/ESC_WT_input_R1.bw \
    output/bigwig_DiffBind_LIB/ESC_WT_input_R2.bw

wiggletools write_bg output/bigwig_DiffBind_LIB/2dN_HET_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/2dN_HET_input_R1.bw \
    output/bigwig_DiffBind_LIB/2dN_HET_input_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/2dN_KO_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/2dN_KO_input_R1.bw \
    output/bigwig_DiffBind_LIB/2dN_KO_input_R2.bw
wiggletools write_bg output/bigwig_DiffBind_LIB/2dN_WT_input_median.bedGraph \
    median output/bigwig_DiffBind_LIB/2dN_WT_input_R1.bw \
    output/bigwig_DiffBind_LIB/2dN_WT_input_R2.bw


# Sort the bedgraph 
bedtools sort -i output/bigwig_DiffBind_LIB/NPC_HET_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/NPC_HET_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/NPC_KO_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/NPC_KO_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/NPC_WT_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/NPC_WT_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/ESC_HET_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/ESC_HET_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/ESC_KO_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/ESC_KO_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/ESC_WT_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/ESC_WT_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/2dN_HET_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/2dN_HET_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/2dN_KO_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/2dN_KO_H3K27me3_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/2dN_WT_H3K27me3_median.bedGraph > \
    output/bigwig_DiffBind_LIB/2dN_WT_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_DiffBind_LIB/NPC_HET_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/NPC_HET_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/NPC_KO_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/NPC_KO_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/NPC_WT_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/NPC_WT_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/ESC_HET_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/ESC_HET_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/ESC_KO_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/ESC_KO_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/ESC_WT_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/ESC_WT_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/2dN_HET_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/2dN_HET_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/2dN_KO_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/2dN_KO_input_median.sorted.bedGraph
bedtools sort -i output/bigwig_DiffBind_LIB/2dN_WT_input_median.bedGraph > \
    output/bigwig_DiffBind_LIB/2dN_WT_input_median.sorted.bedGraph
    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_DiffBind_LIB/NPC_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/NPC_HET_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/NPC_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/NPC_KO_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/NPC_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/NPC_WT_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/ESC_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/ESC_HET_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/ESC_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/ESC_KO_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/ESC_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/ESC_WT_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/2dN_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/2dN_HET_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/2dN_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/2dN_KO_H3K27me3_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/2dN_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/2dN_WT_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_DiffBind_LIB/NPC_HET_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/NPC_HET_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/NPC_KO_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/NPC_KO_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/NPC_WT_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/NPC_WT_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/ESC_HET_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/ESC_HET_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/ESC_KO_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/ESC_KO_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/ESC_WT_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/ESC_WT_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/2dN_HET_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/2dN_HET_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/2dN_KO_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/2dN_KO_input_median.bw
bedGraphToBigWig output/bigwig_DiffBind_LIB/2dN_WT_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_DiffBind_LIB/2dN_WT_input_median.bw