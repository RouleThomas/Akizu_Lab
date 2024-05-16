#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_EZH1cs_peaks.broadPeak \
    -S output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_TMM/PSCWTvsKOEF1aEZH1EZH2TMM-s1-rep0.bw output/THOR/THOR_PSC_WTvsKOEF1aEZH1_EZH2_TMM/PSCWTvsKOEF1aEZH1EZH2TMM-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall1.pdf \
    --samplesLabel "EZH2_WT" "EZH2_KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 10 10 


plotHeatmap -m output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall7.pdf \
    --samplesLabel "EZH2_WT" "EZH2_KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 20 20 




plotProfile -m output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103_profile.pdf \
    --samplesLabel "EZH2_WT" "EZH2_KOEF1aEZH1" \
    --perGroup \
    --colors black blue \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



plotHeatmap -m output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall.pdf \
    --samplesLabel "EZH2_WT" "EZH2_KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall2.pdf \
    --samplesLabel "EZH2_WT" "EZH2_KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall3.pdf \
    --samplesLabel "EZH2_WT" "EZH2_KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4




plotHeatmap -m output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_EZH2_KOEF1aEZH1006vsWT006_THOR_TMM_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall4.pdf \
    --samplesLabel "EZH2_WT" "EZH2_KOEF1aEZH1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142" \
    --colorNumber 10
