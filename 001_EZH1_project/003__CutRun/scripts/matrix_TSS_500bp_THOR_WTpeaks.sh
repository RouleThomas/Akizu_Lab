#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 500 -a 500 \
    -R output/macs2/broad_blacklist_qval2.30103/8wN_WT_H3K27me3_pool_peaks.broadPeak  \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_500bp_THOR_WTpeaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_500bp_THOR_WTpeaks.gz \
    -out output/deeptools/matrix_TSS_500bp_THOR_WTpeaks_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList white,grey,black \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --refPointLabel 0
