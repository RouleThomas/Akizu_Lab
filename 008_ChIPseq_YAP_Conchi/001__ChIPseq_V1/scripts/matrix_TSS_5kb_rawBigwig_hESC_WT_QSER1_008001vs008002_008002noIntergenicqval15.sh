#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../002__ChIPseq_Dixon2021/meta/ENCFF159KBI_macs2_hESC_WT_QSER1FLAG_pool_qval15_noIntergenic.gtf \
    -S output/bigwig/hESC_WT_QSER1_median.bw ../002__ChIPseq_Dixon2021/output/bigwig/hESC_WT_QSER1FLAG_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15_heatmap_colorSmall1.pdf \
    --samplesLabel "QSER1" "QSER1FLAG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 5 80


plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15_heatmap_colorSmall.pdf \
    --samplesLabel "QSER1" "QSER1FLAG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15_heatmap_colorSmall2.pdf \
    --samplesLabel "QSER1" "QSER1FLAG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15_heatmap_colorSmall3.pdf \
    --samplesLabel "QSER1" "QSER1FLAG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15_profile.pdf \
    --samplesLabel "QSER1" "QSER1FLAG" \
    --perGroup \
    --colors black darkblue \
    --refPointLabel "TSS" \
    -T "Read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_QSER1_008001vs008002_008002noIntergenicqval15_heatmap_colorSmall4.pdf \
    --samplesLabel "QSER1" "QSER1FLAG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



