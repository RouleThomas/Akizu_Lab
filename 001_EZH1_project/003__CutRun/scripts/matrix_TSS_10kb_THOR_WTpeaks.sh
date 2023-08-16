#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_pool_peaks.broadPeak  \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_THOR_WTpeaks.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_WTpeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_WTpeaks_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --refPointLabel 0


plotProfile -m output/deeptools/matrix_TSS_10kb_THOR_WTpeaks.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_WTpeaks_profile.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --perGroup \
    --colors black blue red \
    --refPointLabel "0" \
    -T "H3K27me3 read density" \
    -z ""

