#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_NPC_SUZ12/THOR_qval10_positive.bed output/THOR/THOR_NPC_SUZ12/THOR_qval10_negative.bed \
    -S output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s1-rep0.bw output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_THOR_SUZ12_q10_peaks_positive_negative.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12_q10_peaks_positive_negative.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12_q10_peaks_positive_negative_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2


plotProfile -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12_q10_peaks_positive_negative.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12_q10_peaks_positive_negative_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "0" \
    -T "SUZ12 read density" \
    -z ""