#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/sicer2/window200gap600fdr01evalue1000_initialBigwig/NPC_KO_H3K27me3_005008-increased.bed output/sicer2/window200gap600fdr01evalue1000_initialBigwig/NPC_KO_H3K27me3_005008-decreased.bed \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_initialBigwig_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb-SICER2-NPC_WTKO_H3K27me3_window200gap600fdr01evalue1000_initialBigwig-peak.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb-SICER2-NPC_WTKO_H3K27me3_window200gap600fdr01evalue1000_initialBigwig-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-SICER2-NPC_WTKO_H3K27me3_window200gap600fdr01evalue1000_initialBigwig-peak_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 10 10


plotHeatmap -m output/deeptools/matrix_TSS_10kb-SICER2-NPC_WTKO_H3K27me3_window200gap600fdr01evalue1000_initialBigwig-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-SICER2-NPC_WTKO_H3K27me3_window200gap600fdr01evalue1000_initialBigwig-peak_heatmap1.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 5 5

plotHeatmap -m output/deeptools/matrix_TSS_10kb-SICER2-NPC_WTKO_H3K27me3_window200gap600fdr01evalue1000_initialBigwig-peak.gz \
    -out output/deeptools/matrix_TSS_10kb-SICER2-NPC_WTKO_H3K27me3_window200gap600fdr01evalue1000_initialBigwig-peak_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 2 2
