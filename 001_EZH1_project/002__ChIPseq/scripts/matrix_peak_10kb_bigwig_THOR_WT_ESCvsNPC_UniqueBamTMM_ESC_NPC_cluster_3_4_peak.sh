#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4_sort.bed \
    -S output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1_median.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_peak_10kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4_sort.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_peak_10kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4_sort.gz \
    -out output/deeptools/matrix_peak_10kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4_3_4_sort_heatmap.png \
    --colorMap bwr \
    --refPointLabel center

