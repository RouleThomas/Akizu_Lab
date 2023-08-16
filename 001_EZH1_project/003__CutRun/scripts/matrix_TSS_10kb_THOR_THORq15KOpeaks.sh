#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_WTvsKO_unique_Keepdup/THOR_qval15.bed \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kp_THORq15KOpeaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kp_THORq15KOpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15


plotProfile -m output/deeptools/matrix_TSS_10kp_THORq15KOpeaks.gz \
    -out output/deeptools/matrix_TSS_10kp_THORq15KOpeaks_profile.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --perGroup \
    --colors black blue red \
    --refPointLabel "0" \
    -T "H3K27me3 read density" \
    -z ""

