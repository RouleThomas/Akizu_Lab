#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



computeMatrix scale-regions \
    -b 2000 -a 2000 \
    -R ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster4.gtf \
    -S output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_WT_H3K27me3_median.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_WT_H3K27me3_median.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/ESC_HET_H3K27me3_median.bw output/bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF/NPC_HET_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_2kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_HETcluster4.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_gene_2kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_HETcluster4.gz \
    -out output/deeptools/matrix_gene_2kb_bigwig_ChIPseqSpikeInFree_BamToBedToBigwig_uniqueSF_HETcluster4_heatmap.png \
    --colorMap bwr

