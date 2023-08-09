#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



computeMatrix scale-regions \
    -b 2000 -a 2000 \
    -R ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster3.gtf \
    -S output/bigwig_DiffBind_TMM/ESC_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/NPC_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/ESC_HET_H3K27me3_median.bw output/bigwig_DiffBind_TMM/NPC_HET_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_2kb_bigwig_DiffBind_TMM_HETcluster3.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_gene_2kb_bigwig_DiffBind_TMM_HETcluster3.gz \
    -out output/deeptools/matrix_gene_2kb_bigwig_DiffBind_TMM_HETcluster3_heatmap.png \
    --colorMap bwr \
    --refPointLabel center

