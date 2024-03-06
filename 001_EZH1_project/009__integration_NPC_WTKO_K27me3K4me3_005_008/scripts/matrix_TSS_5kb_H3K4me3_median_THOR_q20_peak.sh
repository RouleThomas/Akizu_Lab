#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/THOR/THOR_NPC_WTvsKO_H3K4me3/THOR_qval20.bed \
    -S output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K4me3/NPCWTvsKOH3K4me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_H3K4me3_median_THOR_q20_peak.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb_H3K4me3_median_THOR_q20_peak.gz \
    -out output/deeptools/matrix_TSS_5kb_H3K4me3_median_THOR_q20_peak_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


