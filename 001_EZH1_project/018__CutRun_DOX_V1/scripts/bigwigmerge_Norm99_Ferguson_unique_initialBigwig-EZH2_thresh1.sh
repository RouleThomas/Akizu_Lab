#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bedGraph \
    median output/bigwig_Ferguson/ESC_WT_EZH2_R1_noXchr_unique_norm99_initialBigwig_thresh1.bw \
    output/bigwig_Ferguson/ESC_WT_EZH2_R2_noXchr_unique_norm99_initialBigwig_thresh1.bw \
    output/bigwig_Ferguson/ESC_WT_EZH2_R3_noXchr_unique_norm99_initialBigwig_thresh1.bw 
wiggletools write_bg output/bigwig_Ferguson/ESC_KO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bedGraph \
    median output/bigwig_Ferguson/ESC_KO_EZH2_R1_noXchr_unique_norm99_initialBigwig_thresh1.bw \
    output/bigwig_Ferguson/ESC_KO_EZH2_R2_noXchr_unique_norm99_initialBigwig_thresh1.bw \
    output/bigwig_Ferguson/ESC_KO_EZH2_R3_noXchr_unique_norm99_initialBigwig_thresh1.bw
wiggletools write_bg output/bigwig_Ferguson/ESC_OEKO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bedGraph \
    median output/bigwig_Ferguson/ESC_OEKO_EZH2_R1_noXchr_unique_norm99_initialBigwig_thresh1.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH2_R2_noXchr_unique_norm99_initialBigwig_thresh1.bw \
    output/bigwig_Ferguson/ESC_OEKO_EZH2_R3_noXchr_unique_norm99_initialBigwig_thresh1.bw

# Sort the bedgraph 
bedtools sort -i output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bedGraph > output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/ESC_KO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bedGraph > output/bigwig_Ferguson/ESC_KO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/ESC_OEKO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bedGraph > output/bigwig_Ferguson/ESC_OEKO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_WT_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw

bedGraphToBigWig output/bigwig_Ferguson/ESC_KO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_KO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw

bedGraphToBigWig output/bigwig_Ferguson/ESC_OEKO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_OEKO_EZH2_noXchr_unique_norm99_initialBigwig_thresh1_median.bw


