#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_positive.bed output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_negative.bed \
    -S output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall1.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr

plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList white,blue 


# interactive

plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_10kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall3.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2




