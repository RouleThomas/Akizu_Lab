#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/hESC_WT_QSER1FLAG_median.unique.dupmark.sorted.bedGraph \
    median output/bigwig/hESC_WT_QSER1FLAG_R1.unique.dupmark.sorted.bw \
    output/bigwig/hESC_WT_QSER1FLAG_R2.unique.dupmark.sorted.bw



# Sort the bedgraph 
bedtools sort -i output/bigwig/hESC_WT_QSER1FLAG_median.bedGraph > \
    output/bigwig/hESC_WT_QSER1FLAG_median.sorted.bedGraph

# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/hESC_WT_QSER1FLAG_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/hESC_WT_QSER1FLAG_median.unique.dupmark.sorted.bw




