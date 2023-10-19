#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster1-24_macs2Format.bed output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_heatmap_kmeans25_cluster25.bed \
    -S output/bigwig/PSC_KOEF1aEZH1_SUZ12.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_EZH1cs.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_H3K27me3.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_cl1-24_25.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_cl1-24_25.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_SUZ12_EZH1cs_KOEF1aEZH1_macs2_broadPeak_cl1-24_25_heatmap.pdf \
    --samplesLabel "SUZ12" "EZH1" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'plot, heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --zMin 0 0 0 --zMax 20 10 40

