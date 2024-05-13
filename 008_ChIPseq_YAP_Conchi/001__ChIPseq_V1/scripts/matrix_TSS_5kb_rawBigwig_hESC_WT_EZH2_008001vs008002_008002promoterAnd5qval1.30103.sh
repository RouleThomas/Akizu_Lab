#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R ../002__ChIPseq_Dixon2021/meta/ENCFF159KBI_macs2_hESC_WT_EZH2_R1_qval1.30103_promoterAnd5.gtf \
    -S output/bigwig/hESC_WT_EZH2_median.bw ../002__ChIPseq_Dixon2021/output/bigwig/hESC_WT_EZH2_R1.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103_heatmap_colorSmall1.pdf \
    --samplesLabel "Conchi_EZH2" "Dixon_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 5 40



plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103_heatmap_colorSmall.pdf \
    --samplesLabel "Conchi_EZH2" "Dixon_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103_heatmap_colorSmall2.pdf \
    --samplesLabel "Conchi_EZH2" "Dixon_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103_heatmap_colorSmall3.pdf \
    --samplesLabel "Conchi_EZH2" "Dixon_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103_profile.pdf \
    --samplesLabel "Conchi_EZH2" "Dixon_EZH2" \
    --perGroup \
    --colors black darkblue \
    --refPointLabel "TSS" \
    -T "Read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_rawBigwig_hESC_WT_EZH2_008001vs008002_008002promoterAnd5qval1.30103_heatmap_colorSmall4.pdf \
    --samplesLabel "Conchi_EZH2" "Dixon_EZH2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



