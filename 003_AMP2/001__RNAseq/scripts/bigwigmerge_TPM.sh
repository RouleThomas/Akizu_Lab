#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/CB_Het_median.bedGraph \
    median output/bigwig/CB14_Het_Aligned.sortedByCoord.out.bw \
    output/bigwig/CB42_Het_Aligned.sortedByCoord.out.bw \
    output/bigwig/CB43_Het_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/CB_KO_median.bedGraph \
    median output/bigwig/CB20_KO_Aligned.sortedByCoord.out.bw \
    output/bigwig/CB38_KO_Aligned.sortedByCoord.out.bw \
    output/bigwig/CB41_KO_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/CT_Het_median.bedGraph \
    median output/bigwig/CT14_Het_Aligned.sortedByCoord.out.bw \
    output/bigwig/CT42_Het_Aligned.sortedByCoord.out.bw \
    output/bigwig/CT43_Het_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/CT_KO_median.bedGraph \
    median output/bigwig/CT20_KO_Aligned.sortedByCoord.out.bw \
    output/bigwig/CT38_KO_Aligned.sortedByCoord.out.bw \
    output/bigwig/CT41_KO_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/HP_Het_median.bedGraph \
    median output/bigwig/HP14_Het_Aligned.sortedByCoord.out.bw \
    output/bigwig/HP42_Het_Aligned.sortedByCoord.out.bw \
    output/bigwig/HP43_Het_Aligned.sortedByCoord.out.bw
wiggletools write_bg output/bigwig/HP_KO_median.bedGraph \
    median output/bigwig/HP20_KO_Aligned.sortedByCoord.out.bw \
    output/bigwig/HP38_KO_Aligned.sortedByCoord.out.bw \
    output/bigwig/HP41_KO_Aligned.sortedByCoord.out.bw


# Sort the bedgraph 
bedtools sort -i output/bigwig/CB_Het_median.bedGraph > output/bigwig/CB_Het_median.sorted.bedGraph
bedtools sort -i output/bigwig/CB_KO_median.bedGraph > output/bigwig/CB_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig/CT_Het_median.bedGraph > output/bigwig/CT_Het_median.sorted.bedGraph
bedtools sort -i output/bigwig/CT_KO_median.bedGraph > output/bigwig/CT_KO_median.sorted.bedGraph
bedtools sort -i output/bigwig/HP_Het_median.bedGraph > output/bigwig/HP_Het_median.sorted.bedGraph
bedtools sort -i output/bigwig/HP_KO_median.bedGraph > output/bigwig/HP_KO_median.sorted.bedGraph



# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/CB_Het_median.sorted.bedGraph \
    ../../Master/meta_mice/STAR_mm10/chrNameLength.txt \
    output/bigwig/CB_Het_median.bw
bedGraphToBigWig output/bigwig/CB_KO_median.sorted.bedGraph \
    ../../Master/meta_mice/STAR_mm10/chrNameLength.txt \
    output/bigwig/CB_KO_median.bw
bedGraphToBigWig output/bigwig/CT_Het_median.sorted.bedGraph \
    ../../Master/meta_mice/STAR_mm10/chrNameLength.txt \
    output/bigwig/CT_Het_median.bw
bedGraphToBigWig output/bigwig/CT_KO_median.sorted.bedGraph \
    ../../Master/meta_mice/STAR_mm10/chrNameLength.txt \
    output/bigwig/CT_KO_median.bw
bedGraphToBigWig output/bigwig/HP_Het_median.sorted.bedGraph \
    ../../Master/meta_mice/STAR_mm10/chrNameLength.txt \
    output/bigwig/HP_Het_median.bw
bedGraphToBigWig output/bigwig/HP_KO_median.sorted.bedGraph \
    ../../Master/meta_mice/STAR_mm10/chrNameLength.txt \
    output/bigwig/HP_KO_median.bw

