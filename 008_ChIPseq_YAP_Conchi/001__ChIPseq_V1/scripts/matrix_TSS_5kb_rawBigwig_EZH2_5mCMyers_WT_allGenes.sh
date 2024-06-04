#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI.gtf \
    -S output/bigwig/hESC_WT_EZH2_median.bw output/bigwig/hESC_WT_QSER1_median.bw ../004__WGBS_Dixon2021/output/ENCODE/ENCFF725YJG.bigWig ../004__WGBS_Dixon2021/output/ENCODE/ENCFF040LKO.bigWig \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes_heatmap_colorSmall1.pdf \
    --samplesLabel "EZH2" "QSER1" "5mC_R1" "5mC_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 1.25 2 1 1


plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "EZH2" "QSER1" "5mC_R1" "5mC_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes_heatmap_colorSmall2.pdf \
    --samplesLabel "EZH2" "QSER1" "5mC_R1" "5mC_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes_heatmap_colorSmall3.pdf \
    --samplesLabel "EZH2" "QSER1" "5mC_R1" "5mC_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes_profile.pdf \
    --samplesLabel "EZH2" "QSER1" "5mC_R1" "5mC_R2" \
    --perGroup \
    --colors blue red yellow green purple \
    --refPointLabel "TSS" \
    -T "Read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes_heatmap_colorSmall4.pdf \
    --samplesLabel "EZH2" "QSER1" "5mC_R1" "5mC_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



plotProfile -m output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes.gz \
      --perGroup \
      --plotType heatmap \
      -out output/deeptools/matrix_TSS_5kb_rawBigwig_EZH2_5mCMyers_WT_allGenes_profileHeatmap.pdf \
      --plotHeight 5 \
      --plotWidth 15



