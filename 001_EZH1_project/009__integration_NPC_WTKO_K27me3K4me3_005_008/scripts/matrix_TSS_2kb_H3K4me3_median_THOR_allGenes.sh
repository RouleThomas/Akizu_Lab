#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 2000 -a 2000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_2kb_H3K4me3_median_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_2kb_H3K4me3_median_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_2kb_H3K4me3_median_THOR_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotProfile -m output/deeptools/matrix_TSS_2kb_H3K4me3_median_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_2kb_H3K4me3_median_THOR_allGenes_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K4me3 read density" \
    -z ""



