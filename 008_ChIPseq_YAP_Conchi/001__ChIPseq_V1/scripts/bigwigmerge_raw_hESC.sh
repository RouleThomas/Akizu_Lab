#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/hESC_WT_EZH2_median.bedGraph \
    median output/bigwig/hESC_WT_EZH2_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_EZH2_R2.unique.dupmark.sorted.bw
wiggletools write_bg output/bigwig/hESC_YAPKO_EZH2_median.bedGraph \
    median output/bigwig/hESC_YAPKO_EZH2_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_YAPKO_EZH2_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/hESC_WT_QSER1_median.bedGraph \
    median output/bigwig/hESC_WT_QSER1_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_QSER1_R2.unique.dupmark.sorted.bw
wiggletools write_bg output/bigwig/hESC_YAPKO_QSER1_median.bedGraph \
    median output/bigwig/hESC_YAPKO_QSER1_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_YAPKO_QSER1_R2.unique.dupmark.sorted.bw



# Sort the bedgraph 
bedtools sort -i output/bigwig/hESC_WT_EZH2_median.bedGraph > \
    output/bigwig/hESC_WT_EZH2_median.sorted.bedGraph
bedtools sort -i output/bigwig/hESC_YAPKO_EZH2_median.bedGraph > \
    output/bigwig/hESC_YAPKO_EZH2_median.sorted.bedGraph    

bedtools sort -i output/bigwig/hESC_WT_QSER1_median.bedGraph > \
    output/bigwig/hESC_WT_QSER1_median.sorted.bedGraph
bedtools sort -i output/bigwig/hESC_YAPKO_QSER1_median.bedGraph > \
    output/bigwig/hESC_YAPKO_QSER1_median.sorted.bedGraph 



# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/hESC_WT_EZH2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_EZH2_median.bw
bedGraphToBigWig output/bigwig/hESC_YAPKO_EZH2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_YAPKO_EZH2_median.bw

  bedGraphToBigWig output/bigwig/hESC_WT_QSER1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_QSER1_median.bw
bedGraphToBigWig output/bigwig/hESC_YAPKO_QSER1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_YAPKO_QSER1_median.bw  
    



