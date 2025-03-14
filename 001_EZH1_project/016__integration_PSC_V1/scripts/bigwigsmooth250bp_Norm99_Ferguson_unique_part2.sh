#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


# Convert back bedGraph to bigwig
bedtools sort -i output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth250bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth250bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth250bp.sorted.bedGraph

bedtools sort -i output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth250bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth250bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth250bp.sorted.bedGraph

bedtools sort -i output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth250bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth250bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth250bp.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth250bp.sorted.bedGraph



bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth250bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth250bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth250bp.bw


bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth250bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KO_EZH2_unique_norm99_median_smooth250bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_unique_norm99_median_smooth250bp.bw


bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth250bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_EZH2_unique_norm99_median_smooth250bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth250bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_unique_norm99_median_smooth250bp.bw





