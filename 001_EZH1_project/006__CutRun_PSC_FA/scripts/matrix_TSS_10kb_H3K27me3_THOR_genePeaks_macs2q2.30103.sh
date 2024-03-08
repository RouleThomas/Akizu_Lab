#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI_macs2_H3K27me3_WTKO_qval2.30103_Promoter_5.gtf \
    -S output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-s1-rep0.bw output/THOR/THOR_PSC_WTvsKO_H3K27me3/PSCWTvsKOH3K27me3-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103_heatmap2.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103_heatmap3.pdf \
    --samplesLabel "WT" "KO" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4



plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_THOR_genePeaks_macs2q2.30103_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""



