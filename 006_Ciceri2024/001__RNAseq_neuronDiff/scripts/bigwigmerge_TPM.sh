#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/100dN_WT_median.bedGraph \
    median output/bigwig/100dN_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/100dN_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/100dN_WT_R3_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/75dN_WT_median.bedGraph \
    median output/bigwig/75dN_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/75dN_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/75dN_WT_R3_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/50dN_WT_median.bedGraph \
    median output/bigwig/50dN_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/50dN_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/50dN_WT_R3_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/25dN_WT_median.bedGraph \
    median output/bigwig/25dN_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/25dN_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/25dN_WT_R3_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/NPC_WT_median.bedGraph \
    median output/bigwig/NPC_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/NPC_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/NPC_WT_R3_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/ESC_WT_median.bedGraph \
    median output/bigwig/ESC_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig/ESC_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig/ESC_WT_R3_Aligned.sortedByCoord.out.bw


# Sort the bedgraph 
bedtools sort -i output/bigwig/100dN_WT_median.bedGraph > output/bigwig/100dN_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig/75dN_WT_median.bedGraph > output/bigwig/75dN_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig/50dN_WT_median.bedGraph > output/bigwig/50dN_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig/25dN_WT_median.bedGraph > output/bigwig/25dN_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig/NPC_WT_median.bedGraph > output/bigwig/NPC_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig/ESC_WT_median.bedGraph > output/bigwig/ESC_WT_median.sorted.bedGraph



    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/100dN_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/100dN_WT_median.bw

bedGraphToBigWig output/bigwig/75dN_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/75dN_WT_median.bw

bedGraphToBigWig output/bigwig/50dN_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/50dN_WT_median.bw

bedGraphToBigWig output/bigwig/25dN_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/25dN_WT_median.bw


bedGraphToBigWig output/bigwig/NPC_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/NPC_WT_median.bw

bedGraphToBigWig output/bigwig/ESC_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/ESC_WT_median.bw
