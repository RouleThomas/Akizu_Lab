#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


# Convert a value below 2 as a 0
awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph


awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph


awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 2) $4=0; print}' \
  output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph







# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_WT_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_WT_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_WT_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bw



bedGraphToBigWig output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_KO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_KO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_KO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bw



bedGraphToBigWig output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R1_noXchr_unique_norm99_initialBigwig_thresh2.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R2_noXchr_unique_norm99_initialBigwig_thresh2.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_OEKO_H3K27me3_R3_noXchr_unique_norm99_initialBigwig_thresh2.bw




