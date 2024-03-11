#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap2.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap Greys \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15



plotProfile -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_profile.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --perGroup \
    --colors black blue red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


# 202309 Naiara tasks
plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap2_color.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2

plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap2_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

# 20240310 Naiara plot Slack




plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap3_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap4_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue \
    --colorNumber 3


plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap4_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue,blue,blue \
    --colorNumber 6


plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap5_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap binary

plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap6_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,Gainsboro,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10


plotHeatmap -m output/deeptools/matrix_TSS2_10kp_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS2_10kp_THOR_allGenes_heatmap7_colorSmall.pdf \
    --samplesLabel "WT" "HET" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10

