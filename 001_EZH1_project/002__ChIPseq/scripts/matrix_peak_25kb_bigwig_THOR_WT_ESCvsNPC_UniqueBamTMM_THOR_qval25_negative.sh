#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 25000 -a 25000 \
    -R output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/THOR_qval25_negative.bed \
    -S output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s1_median.bw output/THOR/THOR_WT_ESCvsNPC_UniqueBamTMM/WTESCvsNPCUniqueBamTMM-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_negative.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_negative.bed




plotHeatmap -m output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_negative.gz \
    -out output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_negative_heatmap.png \
    --colorMap bwr \
    --refPointLabel center


