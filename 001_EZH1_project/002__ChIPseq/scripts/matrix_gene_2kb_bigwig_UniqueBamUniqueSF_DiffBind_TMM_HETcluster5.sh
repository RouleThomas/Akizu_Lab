#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



computeMatrix scale-regions \
    -b 2000 -a 2000 \
    -R ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster5.gtf \
    -S output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/ESC_WT_H3K27me3_R1.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/ESC_WT_H3K27me3_R2.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/NPC_WT_H3K27me3_R1.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/NPC_WT_H3K27me3_R2.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/ESC_HET_H3K27me3_R1.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/ESC_HET_H3K27me3_R2.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/NPC_HET_H3K27me3_R1.bw output/bigwig_UniqueBamUniqueSF_DiffBind_TMM/NPC_HET_H3K27me3_R2.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_2kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HETcluster5.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_gene_2kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HETcluster5.gz \
    -out output/deeptools/matrix_gene_2kb_bigwig_UniqueBamUniqueSF_DiffBind_TMM_HETcluster5_heatmap.png \
    --colorMap bwr

