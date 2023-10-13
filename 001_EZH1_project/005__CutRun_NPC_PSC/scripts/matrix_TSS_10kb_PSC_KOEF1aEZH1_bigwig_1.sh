#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig/PSC_KOEF1aEZH1_EZH1cs.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_SUZ12.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_H3K27me3.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_IGG.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_PSC_KOEF1aEZH1_bigwig_1.gz \
    -p 8


plotHeatmap -m output/deeptools/matrix_TSS_10kb_PSC_KOEF1aEZH1_bigwig_1.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_KOEF1aEZH1_bigwig_1_heatmap.png \
    --samplesLabel "EZH1cs" "SUZ12" "H3K27me3" "IGG" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --zMin 0 0 0 0 --zMax 1 2 3 3



plotProfile -m output/deeptools/matrix_TSS_10kb_PSC_KOEF1aEZH1_bigwig_1.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_KOEF1aEZH1_bigwig_1_profile.pdf \
    --samplesLabel "EZH1cs" "SUZ12" "H3K27me3" "IGG" \
    --perGroup \
    --colors darkblue green purple gray \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


