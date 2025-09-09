#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00


# Convert a value below 1 as a 0
awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_WT_EZH1_R1_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_WT_EZH1_R1_unique_norm99_initialBigwig_thresh1.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_WT_EZH1_R2_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_WT_EZH1_R2_unique_norm99_initialBigwig_thresh1.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_WT_EZH1_R3_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_WT_EZH1_R3_unique_norm99_initialBigwig_thresh1.bedGraph


awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_KO_EZH1_R1_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_KO_EZH1_R1_unique_norm99_initialBigwig_thresh1.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_KO_EZH1_R2_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_KO_EZH1_R2_unique_norm99_initialBigwig_thresh1.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_KO_EZH1_R3_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_KO_EZH1_R3_unique_norm99_initialBigwig_thresh1.bedGraph


awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_OEKO_EZH1_R1_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_OEKO_EZH1_R1_unique_norm99_initialBigwig_thresh1.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_OEKO_EZH1_R2_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_OEKO_EZH1_R2_unique_norm99_initialBigwig_thresh1.bedGraph
awk 'BEGIN{OFS="\t"} {if ($4 <= 1) $4=0; print}' \
  output/bigwig_Ferguson/ESC_OEKO_EZH1_R3_unique_norm99_initialBigwig.sorted.bedGraph \
  > output/bigwig_Ferguson/ESC_OEKO_EZH1_R3_unique_norm99_initialBigwig_thresh1.bedGraph







# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_Ferguson/ESC_WT_EZH1_R1_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_WT_EZH1_R1_unique_norm99_initialBigwig_thresh1.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_WT_EZH1_R2_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_WT_EZH1_R2_unique_norm99_initialBigwig_thresh1.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_WT_EZH1_R3_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_WT_EZH1_R3_unique_norm99_initialBigwig_thresh1.bw



bedGraphToBigWig output/bigwig_Ferguson/ESC_KO_EZH1_R1_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_KO_EZH1_R1_unique_norm99_initialBigwig_thresh1.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_KO_EZH1_R2_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_KO_EZH1_R2_unique_norm99_initialBigwig_thresh1.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_KO_EZH1_R3_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_KO_EZH1_R3_unique_norm99_initialBigwig_thresh1.bw



bedGraphToBigWig output/bigwig_Ferguson/ESC_OEKO_EZH1_R1_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R1_unique_norm99_initialBigwig_thresh1.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_OEKO_EZH1_R2_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R2_unique_norm99_initialBigwig_thresh1.bw
bedGraphToBigWig output/bigwig_Ferguson/ESC_OEKO_EZH1_R3_unique_norm99_initialBigwig_thresh1.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_Ferguson/ESC_OEKO_EZH1_R3_unique_norm99_initialBigwig_thresh1.bw




