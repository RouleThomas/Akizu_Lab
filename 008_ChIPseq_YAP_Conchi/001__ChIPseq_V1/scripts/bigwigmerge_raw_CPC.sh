#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/bigwig/CPC_untreated_NR2F2_median.bedGraph \
    median output/bigwig/CPC_untreated_NR2F2_R1.unique.dupmark.sorted.bw \
    output/bigwig/CPC_untreated_NR2F2_R2.unique.dupmark.sorted.bw

wiggletools write_bg output/bigwig/CPC_RA_NR2F2_median.bedGraph \
    median output/bigwig/CPC_RA_NR2F2_R1.unique.dupmark.sorted.bw \
    output/bigwig/CPC_RA_NR2F2_R2.unique.dupmark.sorted.bw



# Sort the bedgraph 
bedtools sort -i output/bigwig/CPC_untreated_NR2F2_median.bedGraph > \
    output/bigwig/CPC_untreated_NR2F2_median.sorted.bedGraph

bedtools sort -i output/bigwig/CPC_RA_NR2F2_median.bedGraph > \
    output/bigwig/CPC_RA_NR2F2_median.sorted.bedGraph

    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig/CPC_untreated_NR2F2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/CPC_untreated_NR2F2_median.bw

bedGraphToBigWig output/bigwig/CPC_RA_NR2F2_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig/CPC_RA_NR2F2_median.bw





