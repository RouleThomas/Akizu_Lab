#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI.gtf \
    -S output/bigwig/hESC_WT_QSER1_median.bw output/bigwig/hESC_WT_EZH2_median.bw ../003__ChIPseq_pluripotency/output/bigwig/hESC_WT_YAP1_R1.unique.dupmark.sorted.bw ../003__ChIPseq_pluripotency/output/bigwig/hESC_WT_TEAD4_R1.unique.dupmark.sorted.bw output/bigwig/hESC_WT_DVL2_R1.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes_heatmap_colorSmall1.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" "TEAD4" "DVL2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 2 1.25 1 1.25 1.25


plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" "TEAD4" "DVL2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes_heatmap_colorSmall2.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" "TEAD4" "DVL2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes_heatmap_colorSmall3.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" "TEAD4" "DVL2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes_profile.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" "TEAD4" "DVL2" \
    --perGroup \
    --colors blue red yellow green purple \
    --refPointLabel "TSS" \
    -T "Read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes_heatmap_colorSmall4.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" "TEAD4" "DVL2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



plotProfile -m output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes.gz \
      --perGroup \
      --plotType heatmap \
      -out output/deeptools/matrix_TSS_5kb_rawBigwig_QSER1EZH2YAP1TEAD4DVL2_WT_allGenes_profileHeatmap.pdf \
      --plotHeight 5 \
      --plotWidth 15



