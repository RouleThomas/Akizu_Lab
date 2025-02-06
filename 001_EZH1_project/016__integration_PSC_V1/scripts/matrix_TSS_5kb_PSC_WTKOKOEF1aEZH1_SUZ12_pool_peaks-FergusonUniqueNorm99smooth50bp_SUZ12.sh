#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks.sorted.merge.bed \
    -S output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median_smooth50bp.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median_smooth50bp.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12_heatmap_colorSmall.pdf \
    --samplesLabel "WT_SUZ12" "KO_SUZ12" "KOEF1aEZH1_SUZ12" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


# interactive
plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12_heatmap_colorSmall1.pdf \
    --samplesLabel "WT_SUZ12" "KO_SUZ12" "KOEF1aEZH1_SUZ12" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 1

plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_SUZ12" "KO_SUZ12" "KOEF1aEZH1_SUZ12" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 0.5


plotProfile -m output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_WTKOKOEF1aEZH1_SUZ12_pool_peaks-FergusonUniqueNorm99smooth50bp_SUZ12_plotProfile1.pdf \
    --samplesLabel "WT_SUZ12" "KO_SUZ12" "KOEF1aEZH1_SUZ12" \
    --colors black red blue \
    --perGroup


