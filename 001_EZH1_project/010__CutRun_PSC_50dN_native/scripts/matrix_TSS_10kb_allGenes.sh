#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig/PSC_WT_H3K27me3.unique.dupmark.sorted.bw output/bigwig/PSC_WT_EZH2.unique.dupmark.sorted.bw output/bigwig/PSC_WT_H3K27me1.unique.dupmark.sorted.bw output/bigwig/PSC_WT_IGG.unique.dupmark.sorted.bw output/bigwig/50dN_WT_H3K27me3.unique.dupmark.sorted.bw output/bigwig/50dN_WT_EZH2.unique.dupmark.sorted.bw output/bigwig/50dN_WT_H3K27me1.unique.dupmark.sorted.bw output/bigwig/50dN_WT_IGG.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_allGenes.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_10kb_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "PSC_H3K27me3" "PSC_EZH2" "PSC_H3K27me1" "PSC_IGG" "50dN_H3K27me3" "50dN_EZH2" "50dN_H3K27me1" "50dN_IGG" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_allGenes_heatmap2.pdf \
    --samplesLabel "PSC_H3K27me3" "PSC_EZH2" "PSC_H3K27me1" "PSC_IGG" "50dN_H3K27me3" "50dN_EZH2" "50dN_H3K27me1" "50dN_IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_allGenes_heatmap3.pdf \
    --samplesLabel "PSC_H3K27me3" "PSC_EZH2" "PSC_H3K27me1" "PSC_IGG" "50dN_H3K27me3" "50dN_EZH2" "50dN_H3K27me1" "50dN_IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotHeatmap -m output/deeptools/matrix_TSS_10kb_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_allGenes_heatmap4.pdf \
    --samplesLabel "PSC_H3K27me3" "PSC_EZH2" "PSC_H3K27me1" "PSC_IGG" "50dN_H3K27me3" "50dN_EZH2" "50dN_H3K27me1" "50dN_IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue,blue,blue,blue,blue,blue,blue \
    --colorNumber 10



plotProfile -m output/deeptools/matrix_TSS_10kb_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_allGenes_profile.pdf \
    --samplesLabel "PSC_H3K27me3" "PSC_EZH2" "PSC_H3K27me1" "PSC_IGG" "50dN_H3K27me3" "50dN_EZH2" "50dN_H3K27me1" "50dN_IGG" \
    --perGroup \
    --colors black red grey white white white white white \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""




plotHeatmap -m output/deeptools/matrix_TSS_10kb_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_allGenes_heatmap5.pdf \
    --samplesLabel "PSC_H3K27me3" "PSC_EZH2" "PSC_H3K27me1" "PSC_IGG" "50dN_H3K27me3" "50dN_EZH2" "50dN_H3K27me1" "50dN_IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10

