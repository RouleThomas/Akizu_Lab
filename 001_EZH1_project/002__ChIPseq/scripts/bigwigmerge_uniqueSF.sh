#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=100:00:00


# Merge bigwig and compute median into a bedgraph

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_HET_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_HET_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_HET_H3K27me3_R2.bw
    
wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_KO_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_KO_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_KO_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_WT_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_WT_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_WT_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_KO_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_KO_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_KO_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_R2.bw

wiggletools write_bg output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_median.bedGraph \
    median output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_R1.bw \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_R2.bw


# Sort the bedgraph 
bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_HET_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_HET_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_KO_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_KO_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_WT_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_WT_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_KO_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_KO_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_median.sorted.bedGraph

bedtools sort -i output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_median.bedGraph > \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_median.sorted.bedGraph

    
# Convert bedgraph to bigwig
bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_HET_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_KO_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/2dN_WT_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_KO_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_KO_H3K27me3_median.bw

bedGraphToBigWig output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_median.sorted.bedGraph \
    ../../Master/meta/GRCh38_chrom_sizes.tab \
    output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_median.bw