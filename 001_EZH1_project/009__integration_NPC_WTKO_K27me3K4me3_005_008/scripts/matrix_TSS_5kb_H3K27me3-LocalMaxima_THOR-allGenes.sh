#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_H3K27me3-LocalMaxima_THOR-allGenes_heatmap.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" "WT_THOR" "KO_THOR" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 10 10 100 100


