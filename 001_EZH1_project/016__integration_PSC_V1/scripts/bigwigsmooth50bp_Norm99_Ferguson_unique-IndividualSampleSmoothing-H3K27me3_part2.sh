#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


# Convert back bedGraph to bigwig
bedtools sort -i output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_006R_unique_norm99_smooth50bp.bw

bedtools sort -i output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_010R_unique_norm99_smooth50bp.bw


bedtools sort -i output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_WT_H3K27me3_013R1_unique_norm99_smooth50bp.bw






bedtools sort -i output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_006R_unique_norm99_smooth50bp.bw

bedtools sort -i output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_013R1_unique_norm99_smooth50bp.bw

bedtools sort -i output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KO_H3K27me3_014R2_unique_norm99_smooth50bp.bw





bedtools sort -i output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_005R_unique_norm99_smooth50bp.bw

bedtools sort -i output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_006R_unique_norm99_smooth50bp.bw

bedtools sort -i output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_smooth50bp.bedGraph > output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_smooth50bp.sorted.bedGraph
bedGraphToBigWig \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_smooth50bp.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/PSC_KOEF1aEZH1_H3K27me3_013R1_unique_norm99_smooth50bp.bw

