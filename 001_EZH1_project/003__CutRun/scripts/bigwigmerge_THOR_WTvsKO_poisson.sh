#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1_median.bedGraph \
    median output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1-rep0.bw \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1-rep1.bw \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1-rep2.bw \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1-rep3.bw
wiggletools write_bg output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2_median.bedGraph \
    median output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2-rep0.bw \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2-rep1.bw \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2-rep2.bw \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2-rep3.bw  


# Sort the bedgraph 
bedtools sort -i output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1_median.bedGraph > \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1_median.sorted.bedGraph
bedtools sort -i output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2_median.bedGraph > \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2_median.sorted.bedGraph


    
# Convert bedgraph to bigwig
bedGraphToBigWig output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s1_median.bw
bedGraphToBigWig output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_WTvsKO_poisson/WTvsKO-s2_median.bw
