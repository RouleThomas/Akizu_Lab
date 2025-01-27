#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00





# Convert to bedGraph and smooth
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bedGraph
multiBigwigSummary bins -b output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median.bw --binSize 50 -o output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.npz --outRawCounts output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.bedGraph

# Convert back bedGraph to bigwig
bedtools sort -i output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.sorted.bedGraph
bedtools sort -i output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.sorted.bedGraph


bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_EZH2_unique_norm99_median_smooth50bp.bw
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_unique_norm99_median_smooth50bp.bw


