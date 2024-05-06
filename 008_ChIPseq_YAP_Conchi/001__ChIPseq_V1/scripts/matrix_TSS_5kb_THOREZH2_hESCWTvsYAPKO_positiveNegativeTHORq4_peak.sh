#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_positive.bed output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_negative.bed \
    -S output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall1.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr


plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall3.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_profile.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --perGroup \
    --colors black darkred \
    --refPointLabel "center" \
    -T "Read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH2_hESCWTvsYAPKO_positiveNegativeTHORq4_peak_heatmap_colorSmall4.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



