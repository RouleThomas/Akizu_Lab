#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 25000 -a 25000 \
    -R output/macs2_unique/broad_blacklist_qval2.30103/ESC_NPC_WT_H3K27me3_pool_peaks_sort.broadPeak \
    -S output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1_median.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak.bed




plotHeatmap -m output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak.gz \
    -out output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap.png \
    --colorMap bwr \
    --refPointLabel center

plotHeatmap -m output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak.gz \
    -out output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4.png \
    --colorMap bwr \
    --refPointLabel center \
    --kmeans 4 \
    --outFileSortedRegions output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_ESC_NPC_peak_heatmap_cluster4.bed

