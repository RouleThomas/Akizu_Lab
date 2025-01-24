#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_Ferguson/PSC_WT_SUZ12_unique_norm99_median.bw output/bigwig_Ferguson/PSC_KO_SUZ12_unique_norm99_median.bw output/bigwig_Ferguson/PSC_KOEF1aEZH1_SUZ12_unique_norm99_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-FergusonUniqueNorm99.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-FergusonUniqueNorm99.gz \
    -out output/deeptools/matrix_TSS_5kb_PSC_SUZ12_WTKOKOEF1aEZH1-FergusonUniqueNorm99_heatmap_colorSmall.pdf \
    --samplesLabel "WT_SUZ12" "KO_SUZ12" "KOEF1aEZH1_SUZ12" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


