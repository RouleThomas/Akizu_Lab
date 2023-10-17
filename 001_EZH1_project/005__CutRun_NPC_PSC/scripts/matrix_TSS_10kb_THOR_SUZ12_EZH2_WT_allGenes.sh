#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s1-rep0.bw output/THOR/THOR_NPC_EZH2/NPCEZH2-s1-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes_heatmap.pdf \
    --samplesLabel "SUZ12" "EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2

plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes_heatmap_kmeans2.pdf \
    --samplesLabel "SUZ12" "EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 2

plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes_heatmap_kmeans4.pdf \
    --samplesLabel "SUZ12" "EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 4

plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_allGenes_heatmap_kmeans8.pdf \
    --samplesLabel "SUZ12" "EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --kmeans 8

