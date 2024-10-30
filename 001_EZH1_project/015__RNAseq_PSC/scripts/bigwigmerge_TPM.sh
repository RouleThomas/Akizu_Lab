#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/PSC_WT_median.bedGraph \
    median output/bigwig/PSC_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_WT_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig/PSC_KO_median.bedGraph \
    median output/bigwig/PSC_KO_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KO_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KO_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig/PSC_KOEF1aEZH1_median.bedGraph \
    median output/bigwig/PSC_KOEF1aEZH1_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KOEF1aEZH1_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/PSC_KOEF1aEZH1_R3_Aligned.sortedByCoord.out.bw



# Sort the bedgraph 
bedtools sort -i output/bigwig/PSC_WT_median.bedGraph > output/bigwig/PSC_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig/PSC_KO_median.bedGraph > output/bigwig/PSC_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig/PSC_KOEF1aEZH1_median.bedGraph > output/bigwig/PSC_KOEF1aEZH1_median.sorted.bedGraph

    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/PSC_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/PSC_WT_median.bw

bedGraphToBigWig output/bigwig/PSC_KO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/PSC_KO_median.bw

bedGraphToBigWig output/bigwig/PSC_KOEF1aEZH1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/PSC_KOEF1aEZH1_median.bw