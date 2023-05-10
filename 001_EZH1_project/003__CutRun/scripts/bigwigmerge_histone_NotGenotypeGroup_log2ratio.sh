#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00



# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_median.bedGraph \
    median output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_R1_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_R2_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_R3_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_R4_log2ratio.bw


wiggletools write_bg output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_median.bedGraph \
    median output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_R1_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_R2_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_R3_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_R4_log2ratio.bw

  
wiggletools write_bg output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_median.bedGraph \
    median output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_R1_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_R2_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_R3_log2ratio.bw \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_R4_log2ratio.bw   





# Sort the bedgraph 
bedtools sort -i output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_median.bedGraph > \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_median.sorted.bedGraph


bedtools sort -i output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_median.bedGraph > \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_median.sorted.bedGraph


bedtools sort -i output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_median.bedGraph > \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_median.sorted.bedGraph




    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_HET_H3K27me3_log2ratio_median.bw


bedGraphToBigWig output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_KO_H3K27me3_log2ratio_median.bw


bedGraphToBigWig output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_histone_NotGenotypeGroup_IggNorm_log2ratio/8wN_WT_H3K27me3_log2ratio_median.bw


