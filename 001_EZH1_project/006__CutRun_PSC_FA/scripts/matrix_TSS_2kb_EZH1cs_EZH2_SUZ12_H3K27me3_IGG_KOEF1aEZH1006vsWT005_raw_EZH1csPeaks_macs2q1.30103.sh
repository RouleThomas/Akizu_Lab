#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 2000 -a 2000 \
    -R output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_EZH2_peaks.broadPeak \
    -S output/bigwig/PSC_KOEF1aEZH1_EZH1cs.unique.dupmark.sorted.bw output/bigwig/PSC_WT_EZH1cs.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_EZH2.unique.dupmark.sorted.bw output/bigwig/PSC_WT_EZH2.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_SUZ12.unique.dupmark.sorted.bw output/bigwig/PSC_WT_SUZ12.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw output/bigwig/PSC_WT_H3K27me3.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_IGG.unique.dupmark.sorted.bw output/bigwig/PSC_WT_IGG.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.gz \
    -p 6


plotHeatmap -m output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall1.pdf \
    --samplesLabel "EZH1cs_KOEF1aEZH1" "EZH1cs_WT" "EZH2_KOEF1aEZH1" "EZH2_WT" "SUZ12_KOEF1aEZH1" "SUZ12_WT" "H3K27me3_KOEF1aEZH1" "H3K27me3_WT" "IGG_KOEF1aEZH1" "IGG_WT" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr


plotHeatmap -m output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall.pdf \
    --samplesLabel "EZH1cs_KOEF1aEZH1" "EZH1cs_WT" "EZH2_KOEF1aEZH1" "EZH2_WT" "SUZ12_KOEF1aEZH1" "SUZ12_WT" "H3K27me3_KOEF1aEZH1" "H3K27me3_WT" "IGG_KOEF1aEZH1" "IGG_WT" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall2.pdf \
    --samplesLabel "EZH1cs_KOEF1aEZH1" "EZH1cs_WT" "EZH2_KOEF1aEZH1" "EZH2_WT" "SUZ12_KOEF1aEZH1" "SUZ12_WT" "H3K27me3_KOEF1aEZH1" "H3K27me3_WT" "IGG_KOEF1aEZH1" "IGG_WT" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall3.pdf \
    --samplesLabel "EZH1cs_KOEF1aEZH1" "EZH1cs_WT" "EZH2_KOEF1aEZH1" "EZH2_WT" "SUZ12_KOEF1aEZH1" "SUZ12_WT" "H3K27me3_KOEF1aEZH1" "H3K27me3_WT" "IGG_KOEF1aEZH1" "IGG_WT" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103_profile.pdf \
    --samplesLabel "EZH1cs_KOEF1aEZH1" "EZH1cs_WT" "EZH2_KOEF1aEZH1" "EZH2_WT" "SUZ12_KOEF1aEZH1" "SUZ12_WT" "H3K27me3_KOEF1aEZH1" "H3K27me3_WT" "IGG_KOEF1aEZH1" "IGG_WT" \
    --perGroup \
    --colors green purple blue red grey \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_2kb_EZH1cs_EZH2_SUZ12_H3K27me3_IGG_KOEF1aEZH1006vsWT005_raw_EZH1csPeaks_macs2q1.30103_heatmap_colorSmall4.pdf \
    --samplesLabel "EZH1cs_KOEF1aEZH1" "EZH1cs_WT" "EZH2_KOEF1aEZH1" "EZH2_WT" "SUZ12_KOEF1aEZH1" "SUZ12_WT" "H3K27me3_KOEF1aEZH1" "H3K27me3_WT" "IGG_KOEF1aEZH1" "IGG_WT" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10
