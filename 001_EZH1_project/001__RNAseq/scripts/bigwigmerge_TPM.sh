#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_hg38/ESC_WT_median.bedGraph \
    median output/bigwig_hg38/ESC_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/ESC_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/ESC_WT_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/ESC_KO_median.bedGraph \
    median output/bigwig_hg38/ESC_KO_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/ESC_KO_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/ESC_KO_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/ESC_HET_median.bedGraph \
    median output/bigwig_hg38/ESC_HET_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/ESC_HET_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/ESC_HET_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/NPC_WT_median.bedGraph \
    median output/bigwig_hg38/NPC_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/NPC_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/NPC_WT_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/NPC_KO_median.bedGraph \
    median output/bigwig_hg38/NPC_KO_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/NPC_KO_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/NPC_KO_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/NPC_HET_median.bedGraph \
    median output/bigwig_hg38/NPC_HET_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/NPC_HET_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/NPC_HET_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/2dN_WT_median.bedGraph \
    median output/bigwig_hg38/2dN_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/2dN_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/2dN_WT_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/2dN_KO_median.bedGraph \
    median output/bigwig_hg38/2dN_KO_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/2dN_KO_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/2dN_KO_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/2dN_HET_median.bedGraph \
    median output/bigwig_hg38/2dN_HET_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/2dN_HET_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/2dN_HET_R3_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/4wN_WT_median.bedGraph \
    median output/bigwig_hg38/4wN_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/4wN_WT_R2_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/4wN_KO_median.bedGraph \
    median output/bigwig_hg38/4wN_KO_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/4wN_KO_R2_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/4wN_HET_median.bedGraph \
    median output/bigwig_hg38/4wN_HET_R3_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/4wN_HET_R4_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/4wN_iPSCWT_median.bedGraph \
    median output/bigwig_hg38/4wN_iPSCWT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/4wN_iPSCWT_R2_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/4wN_iPSCpatient_median.bedGraph \
    median output/bigwig_hg38/4wN_iPSCpatient_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/4wN_iPSCpatient_R2_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/8wN_WT_median.bedGraph \
    median output/bigwig_hg38/8wN_WT_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_WT_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_WT_R3_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_WT_R4_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/8wN_KO_median.bedGraph \
    median output/bigwig_hg38/8wN_KO_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_KO_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_KO_R3_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_KO_R4_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/8wN_HET_median.bedGraph \
    median output/bigwig_hg38/8wN_HET_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_HET_R2_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_HET_R3_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_HET_R4_Aligned.sortedByCoord.out.bw

wiggletools write_bg output/bigwig_hg38/8wN_iPSCpatient_median.bedGraph \
    median output/bigwig_hg38/8wN_iPSCpatient_R1_Aligned.sortedByCoord.out.bw \
    output/bigwig_hg38/8wN_iPSCpatient_R2_Aligned.sortedByCoord.out.bw
    
    


# Sort the bedgraph 
bedtools sort -i output/bigwig_hg38/ESC_WT_median.bedGraph > output/bigwig_hg38/ESC_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/ESC_KO_median.bedGraph > output/bigwig_hg38/ESC_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/ESC_HET_median.bedGraph > output/bigwig_hg38/ESC_HET_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/NPC_WT_median.bedGraph > output/bigwig_hg38/NPC_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/NPC_KO_median.bedGraph > output/bigwig_hg38/NPC_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/NPC_HET_median.bedGraph > output/bigwig_hg38/NPC_HET_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/2dN_WT_median.bedGraph > output/bigwig_hg38/2dN_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/2dN_KO_median.bedGraph > output/bigwig_hg38/2dN_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/2dN_HET_median.bedGraph > output/bigwig_hg38/2dN_HET_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/4wN_WT_median.bedGraph > output/bigwig_hg38/4wN_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/4wN_KO_median.bedGraph > output/bigwig_hg38/4wN_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/4wN_HET_median.bedGraph > output/bigwig_hg38/4wN_HET_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/4wN_iPSCWT_median.bedGraph > output/bigwig_hg38/4wN_iPSCWT_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/4wN_iPSCpatient_median.bedGraph > output/bigwig_hg38/4wN_iPSCpatient_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/8wN_WT_median.bedGraph > output/bigwig_hg38/8wN_WT_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/8wN_KO_median.bedGraph > output/bigwig_hg38/8wN_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/8wN_HET_median.bedGraph > output/bigwig_hg38/8wN_HET_median.sorted.bedGraph
bedtools sort -i output/bigwig_hg38/8wN_iPSCpatient_median.bedGraph > output/bigwig_hg38/8wN_iPSCpatient_median.sorted.bedGraph


    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_hg38/ESC_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/ESC_WT_median.bw

bedGraphToBigWig output/bigwig_hg38/ESC_KO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/ESC_KO_median.bw

bedGraphToBigWig output/bigwig_hg38/ESC_HET_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/ESC_HET_median.bw

bedGraphToBigWig output/bigwig_hg38/NPC_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/NPC_WT_median.bw

bedGraphToBigWig output/bigwig_hg38/NPC_KO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/NPC_KO_median.bw

bedGraphToBigWig output/bigwig_hg38/NPC_HET_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/NPC_HET_median.bw

bedGraphToBigWig output/bigwig_hg38/2dN_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/2dN_WT_median.bw

bedGraphToBigWig output/bigwig_hg38/2dN_KO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/2dN_KO_median.bw

bedGraphToBigWig output/bigwig_hg38/2dN_HET_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/2dN_HET_median.bw

bedGraphToBigWig output/bigwig_hg38/4wN_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/4wN_WT_median.bw

bedGraphToBigWig output/bigwig_hg38/4wN_KO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/4wN_KO_median.bw

bedGraphToBigWig output/bigwig_hg38/4wN_HET_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/4wN_HET_median.bw

bedGraphToBigWig output/bigwig_hg38/4wN_iPSCWT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/4wN_iPSCWT_median.bw

bedGraphToBigWig output/bigwig_hg38/4wN_iPSCpatient_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/4wN_iPSCpatient_median.bw
    
bedGraphToBigWig output/bigwig_hg38/8wN_WT_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/8wN_WT_median.bw

bedGraphToBigWig output/bigwig_hg38/8wN_KO_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/8wN_KO_median.bw

bedGraphToBigWig output/bigwig_hg38/8wN_HET_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/8wN_HET_median.bw

bedGraphToBigWig output/bigwig_hg38/8wN_iPSCpatient_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_hg38/8wN_iPSCpatient_median.bw

