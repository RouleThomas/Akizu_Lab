#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_peak.gtf \
    -S output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_HET_H3K27me3_median.bw output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks.gz \
    -p 6 \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks.bed



plotProfile -m output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks.gz \
    -out output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_profile.png \
    --perGroup \
    --colors black blue red \
    --samplesLabel "WT" "HET" "KO" \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""
plotHeatmap -m output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks.gz \
    -out output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_heatmap_kmeans6.png \
    --perGroup \
    --kmeans 6 \
    --samplesLabel "WT" "HET" "KO" \
    --outFileSortedRegions output/deeptools/matrix_TSS_5kb_DiffBind_TMM_peaks_heatmap_kmeans6.txt

