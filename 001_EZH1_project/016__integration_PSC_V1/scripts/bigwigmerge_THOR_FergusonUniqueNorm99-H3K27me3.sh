#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00





# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1_median.bedGraph \
    median output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1-rep0.bw \
    output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1-rep1.bw \
    output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1-rep2.bw
wiggletools write_bg output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2_median.bedGraph \
    median output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2-rep0.bw \
    output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2-rep1.bw \
    output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2-rep2.bw

wiggletools write_bg output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1_median.bedGraph \
    median output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1-rep0.bw \
    output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1-rep1.bw \
    output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1-rep2.bw
wiggletools write_bg output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2_median.bedGraph \
    median output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2-rep0.bw \
    output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2-rep1.bw \
    output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2-rep2.bw
    

# Sort the bedgraph 
bedtools sort -i output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1_median.bedGraph > output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1_median.sorted.bedGraph
bedtools sort -i output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2_median.bedGraph > output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2_median.sorted.bedGraph
bedtools sort -i output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1_median.bedGraph > output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1_median.sorted.bedGraph
bedtools sort -i output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2_median.bedGraph > output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2_median.sorted.bedGraph


# Convert bedgraph to bigwig
bedGraphToBigWig output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s1_median.bw

bedGraphToBigWig output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_PSC_WTvsKO_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOH3K27me3FergusonUniqueNorm99noInput-s2_median.bw

bedGraphToBigWig output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s1_median.bw

bedGraphToBigWig output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_PSC_WTvsKOEF1aEZH1_H3K27me3_FergusonUniqueNorm99_noInput/PSCWTvsKOEF1aEZH1H3K27me3FergusonUniqueNorm99noInput-s2_median.bw
