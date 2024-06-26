#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30_gain.bed output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30_lost.bed \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw ../005__CutRun_NPC_PSC/output/THOR/THOR_NPC_EZH2_TMM/NPCEZH2TMM-s1-rep0.bw ../005__CutRun_NPC_PSC/output/THOR/THOR_NPC_EZH2_TMM/NPCEZH2TMM-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_heatmap.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 120 120 10 10 

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_heatmap7.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 120 120 15 15 

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_heatmap8.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 120 120 7 7 

plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_heatmap1.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2 \
    --zMax 120 120 10 10 




plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_heatmap2.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4 \
    --zMax 120 120 10 10 



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_heatmap3.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10 \
    --zMax 120 120 10 10 



plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_profile.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --colors black darkgrey \
    -T "Read density" \
    --yMax 120 120 10 10 \
    --numPlotsPerRow 2

plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_profile1.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --colors black darkgrey red darkred \
    -T "Read density" \
    --perGroup \
    --yMax 120 120 10 10 \
    --numPlotsPerRow 2

plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_EZH2_median_THORTMM_gainLost_profile2.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" \
    --colors white white black red \
    -T "Read density" \
    --perGroup \
    --yMax 4.5 4.5 \
    --yMin 1 1 \
    --numPlotsPerRow 2 \
    --legendLocation upper-left

