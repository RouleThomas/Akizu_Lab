#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 50000 -a 50000 \
    -R output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30_gain.bed output/THOR/THOR_NPC_WTvsKO_H3K27me3/THOR_qval30_lost.bed \
    -S output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s1_median.bw output/THOR/THOR_NPC_WTvsKO_H3K27me3/NPCWTvsKOH3K27me3-s2_median.bw ../005__CutRun_NPC_PSC/output/THOR/THOR_NPC_EZH2_TMM/NPCEZH2TMM-s1-rep0.bw ../005__CutRun_NPC_PSC/output/THOR/THOR_NPC_EZH2_TMM/NPCEZH2TMM-s2-rep0.bw ../005__CutRun_NPC_PSC/output/THOR/THOR_NPC_SUZ12_TMM/NPCSUZ12TMM-s1-rep0.bw ../005__CutRun_NPC_PSC/output/THOR/THOR_NPC_SUZ12_TMM/NPCSUZ12TMM-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost_heatmap.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" "SUZ12_WT" "SUZ12_KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 120 120 10 10 10 10


plotHeatmap -m output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost_heatmap1.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" "SUZ12_WT" "SUZ12_KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2 \
    --zMax 120 120 10 10 10 10




plotHeatmap -m output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost_heatmap2.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" "SUZ12_WT" "SUZ12_KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4 \
    --zMax 120 120 10 10 10 10



plotHeatmap -m output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost_heatmap3.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" "SUZ12_WT" "SUZ12_KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10 \
    --zMax 120 120 10 10 10 10



plotProfile -m output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost.gz \
    -out output/deeptools/matrix_TSS_50kb_H3K27me3_EZH2_SUZ12_median_THORTMM_gainLost_profile.pdf \
    --samplesLabel "H3K27me3_WT" "H3K27me3_KO" "EZH2_WT" "EZH2_KO" "SUZ12_WT" "SUZ12_KO" \
    --colors black darkgrey \
    -T "Read density" \
    --yMax 120 120 6 6 6 6 \
    --numPlotsPerRow 2



