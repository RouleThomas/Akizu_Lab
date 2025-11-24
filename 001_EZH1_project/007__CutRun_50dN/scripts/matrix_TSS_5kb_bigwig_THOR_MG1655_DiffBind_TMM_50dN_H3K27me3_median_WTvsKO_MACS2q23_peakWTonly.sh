#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval2.30103/50dN_WTQ731E_H3K27me3_pool_peaks.broadPeak \
    -S output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bw \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly.gz \
    -p 8


plotHeatmap -m output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly.gz \
    -out output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly_heatmap.png \
    --samplesLabel "WTQ731E" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotHeatmap -m output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly.gz \
    -out output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly_heatmap1.pdf \
    --samplesLabel "WTQ731E" "KO" \
    --colorMap Blues \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2







plotProfile -m output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly.gz \
    -out output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly_profile.pdf \
    --samplesLabel "WTQ731E" "KO" \
    --perGroup \
    --colors black red \
    -T "H3K27me3 read density" \
    -z "" \
    --plotWidth 4




plotHeatmap -m output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly.gz \
    -out output/deeptools/matrix_TSS_5kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_WTvsKO_MACS2q23_peakWTonly_heatmap3.pdf \
    --samplesLabel "WTQ731E" "KO" \
    --colorMap Blues \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3





