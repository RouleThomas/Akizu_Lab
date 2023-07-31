#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph
wiggletools write_bg output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1.bedGraph \
    median output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1-rep0.bw \
    output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1-rep1.bw

wiggletools write_bg output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2.bedGraph \
    median output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2-rep0.bw \
    output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2-rep1.bw

# Sort the bedgraph 
bedtools sort -i output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1.bedGraph > \
    output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1.sorted.bedGraph

bedtools sort -i output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2.bedGraph > \
    output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2.sorted.bedGraph

# Convert bedgraph to bigwig
bedGraphToBigWig output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1_median.bw

bedGraphToBigWig output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2_median.bw

