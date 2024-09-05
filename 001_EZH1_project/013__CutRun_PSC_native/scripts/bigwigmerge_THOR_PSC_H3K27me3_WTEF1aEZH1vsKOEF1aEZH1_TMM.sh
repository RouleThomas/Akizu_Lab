#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
## WTEF1aEZH1
wiggletools write_bg output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s1_median.bedGraph \
    median output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s1-rep0.bw \
    output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s1-rep1.bw
## KOEF1aEZH1
wiggletools write_bg output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s2_median.bedGraph \
    median output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s2-rep0.bw \
    output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s2-rep1.bw
    




# Sort the bedgraph 
bedtools sort -i output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s1_median.bedGraph > \
    output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s1_median.sorted.bedGraph
bedtools sort -i output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s2_median.bedGraph > \
    output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s2_median.sorted.bedGraph
    

    
# Convert bedgraph to bigwig
bedGraphToBigWig output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s1_median.bw
bedGraphToBigWig output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_PSC_H3K27me3_WTEF1aEZH1vsKOEF1aEZH1/PSCH3K27me3WTEF1aEZH1vsKOEF1aEZH1-s2_median.bw
    
