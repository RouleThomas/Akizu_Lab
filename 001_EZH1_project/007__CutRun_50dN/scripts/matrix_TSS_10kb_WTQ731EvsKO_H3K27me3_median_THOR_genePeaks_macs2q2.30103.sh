#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI_macs2_H3K27me3_WTQ731EKO_pool_qval2.30103_Promoter_5.gtf \
    -S output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103_heatmap3.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4




plotProfile -m output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_WTQ731EvsKO_H3K27me3_median_THOR_genePeaks_macs2q2.30103_heatmap4.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



    