#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI__NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Gain_annot_promoterAnd5.gtf meta/ENCFF159KBI__NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_log2FC1_Lost_annot_promoterAnd5.gtf \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb-DIFFREPS-NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_fc1-gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb-DIFFREPS-NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_fc1-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-DIFFREPS-NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_fc1-gene_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 10 10


plotHeatmap -m output/deeptools/matrix_TSS_5kb-DIFFREPS-NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_fc1-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-DIFFREPS-NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_fc1-gene_heatmap1.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 5 5

plotHeatmap -m output/deeptools/matrix_TSS_5kb-DIFFREPS-NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_fc1-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-DIFFREPS-NPC_WTKO_H3K27me3_merged_intervals_5kb2kb1kb500bp250bp__padj001_nb_pval0001_fc1-gene_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 2 2
