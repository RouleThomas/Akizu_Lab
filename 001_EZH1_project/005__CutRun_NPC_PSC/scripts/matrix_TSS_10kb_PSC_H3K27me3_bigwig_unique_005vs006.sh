#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_unique/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw ../006__CutRun_PSC_FA/output/bigwig/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_bigwig_unique_005vs006.gz \
    -p 8


plotHeatmap -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_bigwig_unique_005vs006.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_bigwig_unique_005vs006_heatmap.png \
    --samplesLabel "native" "FA" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 4 \
    --zMin 0 0 --zMax 3 3



plotProfile -m output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_bigwig_unique_005vs006.gz \
    -out output/deeptools/matrix_TSS_10kb_PSC_H3K27me3_bigwig_unique_005vs006_profile.pdf \
    --samplesLabel "native" "FA" \
    --perGroup \
    --colors gray black \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


