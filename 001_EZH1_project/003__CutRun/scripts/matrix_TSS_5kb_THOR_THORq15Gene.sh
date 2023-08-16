#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_THOR_WTvsHETKO_unique_Keepdup_qval15_Promoter_5.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kp_THORq15Gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kp_THORq15Gene.gz \
    -out output/deeptools/matrix_TSS_5kp_THORq15Gene_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15
