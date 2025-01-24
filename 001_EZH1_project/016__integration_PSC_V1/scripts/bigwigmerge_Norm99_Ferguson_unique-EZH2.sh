#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00




# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median.bedGraph \
    median output/bigwig_Ferguson/PSC_WT_EZH2_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_EZH2_010R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_WT_EZH2_014R1_unique_norm99.bw 
wiggletools write_bg output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median.bedGraph \
    median output/bigwig_Ferguson/PSC_KO_EZH2_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_EZH2_014R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KO_EZH2_014R2_unique_norm99.bw
wiggletools write_bg output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median.bedGraph \
    median output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_006R_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_013R1_unique_norm99.bw \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_014R1_unique_norm99.bw

# Sort the bedgraph 
bedtools sort -i output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median.bedGraph > output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median.bedGraph > output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median.bw

bedGraphToBigWig output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median.bw

bedGraphToBigWig output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median.bw


