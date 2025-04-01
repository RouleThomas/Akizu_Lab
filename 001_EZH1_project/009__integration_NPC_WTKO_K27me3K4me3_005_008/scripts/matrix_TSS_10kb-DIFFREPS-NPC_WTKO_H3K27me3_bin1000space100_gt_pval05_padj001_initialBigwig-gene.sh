#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI__NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig_gain_annot_promoterAnd5.gtf meta/ENCFF159KBI__NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig_lost_annot_promoterAnd5.gtf \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb-DIFFREPS-NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig-gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig-gene.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig-gene_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 10 10


plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig-gene.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig-gene_heatmap1.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 5 5

plotHeatmap -m output/deeptools/matrix_TSS_10kb-DIFFREPS-NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig-gene.gz \
    -out output/deeptools/matrix_TSS_10kb-DIFFREPS-NPC_WTKO_H3K27me3_bin1000space100_gt_pval05_padj001_initialBigwig-gene_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 2 2
