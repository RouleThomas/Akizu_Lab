#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00

computeMatrix reference-point --referencePoint TSS \
    -b 1000 -a 1000 \
    -R output/ChIPseeker/THOR_qval15_HET_Gain_DEG_Up.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup.gz \
    -out output/deeptools/matrix_TSS_1kb_bigwig_THOR_WTvsHET_unique_Keepdup_heatmap.png \
    --colorMap bwr

