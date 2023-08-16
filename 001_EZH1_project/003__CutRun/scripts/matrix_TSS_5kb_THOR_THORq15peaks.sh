#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R output/THOR/THOR_WTvsHET_unique_Keepdup/THOR_HETKO_qval15.bed \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kp_THORq15peaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kp_THORq15peaks.gz \
    -out output/deeptools/matrix_TSS_5kp_THORq15peaks_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15
