#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



computeMatrix scale-regions \
    -b 2000 -a 2000 \
    -R ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster3.gtf \
    -S output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s1-rep0.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s1-rep1.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s2-rep0.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s2-rep1.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster3.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster3.gz \
    -out output/deeptools/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster3_heatmap.png \
    --colorMap bwr


# draft



computeMatrix scale-regions \
    -b 2000 -a 2000 \
    -R ../001__RNAseq/output/deseq2_hg38/cluster_gene_rlog_5cl_ESC_NPC_WT_HET_H3K27me3postClustering_cluster5.gtf \
    -S output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s1-rep0.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s1-rep1.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s2-rep0.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamDiffBindTMM/WTESCvsNPCUniqueBamDiffBindTMM-s2-rep1.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster5.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster5.gz \
    -out output/deeptools/matrix_gene_2kb_bigwig_THOR_ESCvsNPC_UniqueBamDiffBindTMM_HETcluster5_heatmap.png \
    --colorMap bwr

