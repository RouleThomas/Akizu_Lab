#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_positive.bed output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/THOR_qval10_negative.bed \
    -S output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq10_positive_negative.gz \
    -p 8


plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq10_positive_negative.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq10_positive_negative_heatmap.png \
    --samplesLabel "WTQ731E" "KOEF1aEZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

plotProfile -m output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq10_positive_negative.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKOEF1aEZH1_THORq10_positive_negative_profile.pdf \
    --samplesLabel "WTQ731E" "KOEF1aEZH1" \
    --perGroup \
    --colors black red \
    --refPointLabel "0" \
    -T "H3K27me3 read density" \
    -z ""


