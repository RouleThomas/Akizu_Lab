#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/hESC_WT_DNMT3A_median.unique.dupmark.sorted.bedGraph \
    median output/bigwig/hESC_WT_DNMT3A_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_DNMT3A_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/hESC_WT_DNMT3B_median.unique.dupmark.sorted.bedGraph \
    median output/bigwig/hESC_WT_DNMT3B_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_DNMT3B_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/hESC_WT_H3K27me3_median.unique.dupmark.sorted.bedGraph \
    median output/bigwig/hESC_WT_H3K27me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_H3K27me3_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/hESC_WT_H3K4me3_median.unique.dupmark.sorted.bedGraph \
    median output/bigwig/hESC_WT_H3K4me3_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_H3K4me3_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/hESC_WT_input_median.unique.dupmark.sorted.bedGraph \
    median output/bigwig/hESC_WT_input_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_input_R2.unique.dupmark.sorted.bw


# Sort the bedgraph 
bedtools sort -i output/bigwig/hESC_WT_DNMT3A_median.bedGraph > \
    output/bigwig/hESC_WT_DNMT3A_median.sorted.bedGraph
bedtools sort -i output/bigwig/hESC_WT_DNMT3B_median.bedGraph > \
    output/bigwig/hESC_WT_DNMT3B_median.sorted.bedGraph

bedtools sort -i output/bigwig/hESC_WT_H3K27me3_median.bedGraph > \
    output/bigwig/hESC_WT_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig/hESC_WT_H3K4me3_median.bedGraph > \
    output/bigwig/hESC_WT_H3K4me3_median.sorted.bedGraph

bedtools sort -i output/bigwig/hESC_WT_input_median.bedGraph > \
    output/bigwig/hESC_WT_input_median.sorted.bedGraph

# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/hESC_WT_DNMT3A_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_DNMT3A_median.unique.dupmark.sorted.bw
bedGraphToBigWig output/bigwig/hESC_WT_DNMT3B_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_DNMT3B_median.unique.dupmark.sorted.bw
bedGraphToBigWig output/bigwig/hESC_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_H3K27me3_median.unique.dupmark.sorted.bw

bedGraphToBigWig output/bigwig/hESC_WT_H3K4me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_H3K4me3_median.unique.dupmark.sorted.bw


bedGraphToBigWig output/bigwig/hESC_WT_input_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_input_median.unique.dupmark.sorted.bw






