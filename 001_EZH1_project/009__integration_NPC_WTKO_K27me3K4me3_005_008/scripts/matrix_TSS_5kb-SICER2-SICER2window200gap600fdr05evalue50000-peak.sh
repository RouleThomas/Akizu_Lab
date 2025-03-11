#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/sicer2/window200gap600fdr05evalue50000/NPC_KO_H3K27me3_005008-increased.bed output/sicer2/window200gap600fdr05evalue50000/NPC_KO_H3K27me3_005008-decreased.bed \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-peak.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-peak.gz \
    -out output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-peak_heatmap.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 10 10




plotHeatmap -m output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-peak.gz \
    -out output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-peak_heatmap1.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 5 5



plotHeatmap -m output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-peak.gz \
    -out output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-peak_heatmap2.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 2 2