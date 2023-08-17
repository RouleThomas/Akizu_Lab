#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_WTvsKO_unique_Keepdup/THOR_qval15_positive.bed output/THOR/THOR_WTvsKO_unique_Keepdup/THOR_qval15_negative.bed \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_positive_negative.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_positive_negative.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_positive_negative_heatmap.png \
    --samplesLabel "WT" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --perGroup \
    --heatmapHeight 15

plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_positive_negative.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_positive_negative_heatmap2.png \
    --samplesLabel "WT" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15

plotProfile -m output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_positive_negative.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_positive_negative_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "0" \
    -T "H3K27me3 read density" \
    -z ""

