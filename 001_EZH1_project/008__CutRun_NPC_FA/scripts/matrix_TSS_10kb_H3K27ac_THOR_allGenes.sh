#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27ac/NPCWTvsKOH3K27ac-s1-rep0.bw output/THOR/THOR_NPC_WTvsKO_H3K27ac/NPCWTvsKOH3K27ac-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27ac_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27ac_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27ac_THOR_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27ac_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27ac_THOR_allGenes_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27ac read density" \
    -z ""



