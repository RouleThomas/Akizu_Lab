#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    --minThreshold 5 \
    -R meta/ENCFF159KBI_peak.gtf \
    -S output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_min5.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_min5.bed



plotProfile -m output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_min5.gz \
    -out output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_min5_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""
plotHeatmap -m output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_min5.gz \
    -out output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_min5_heatmap_kmeans6.png \
    --perGroup \
    --kmeans 6 \
    --colorMap bwr \
    --samplesLabel "WT" "HET" "KO" \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_min5_heatmap_kmeans6.txt

