#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/ESC_WT_median.bedGraph \
    median output/bigwig/ESC_WT_R1.bw \
    output/bigwig/ESC_WT_R2.bw \
    output/bigwig/ESC_WT_R3.bw 
wiggletools write_bg output/bigwig/ESC_KO_median.bedGraph \
    median output/bigwig/ESC_KO_R1.bw \
    output/bigwig/ESC_KO_R2.bw \
    output/bigwig/ESC_KO_R3.bw
wiggletools write_bg output/bigwig/ESC_OEKO_median.bedGraph \
    median output/bigwig/ESC_OEKO_R1.bw \
    output/bigwig/ESC_OEKO_R2.bw \
    output/bigwig/ESC_OEKO_R3.bw

# Sort the bedgraph 
bedtools sort -i output/bigwig/ESC_WT_median.bedGraph > output/bigwig/ESC_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig/ESC_KO_median.bedGraph > output/bigwig/ESC_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig/ESC_OEKO_median.bedGraph > output/bigwig/ESC_OEKO_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/ESC_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/ESC_WT_median.bw

bedGraphToBigWig output/bigwig/ESC_KO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/ESC_KO_median.bw

bedGraphToBigWig output/bigwig/ESC_OEKO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/ESC_OEKO_median.bw


