#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad_blacklist_qval1.30103/PSC_KOEF1aEZH1_SUZ12_peaks.broadPeak \
    -S output/bigwig/PSC_KOEF1aEZH1_SUZ12.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_EZH1cs.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_EZH2.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw output/bigwig/PSC_KOEF1aEZH1_IGG.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103_heatmap_colorSmall1.pdf \
    --samplesLabel "SUZ12" "EZH1cs" "EZH2" "H3K27me3" "IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 10 10 10 80 80


plotHeatmap -m output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103_heatmap_colorSmall.pdf \
    --samplesLabel "SUZ12" "EZH1cs" "EZH2" "H3K27me3" "IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2 \
    --zMax 10 10 10 60 60

plotHeatmap -m output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103_heatmap_colorSmall2.pdf \
    --samplesLabel "SUZ12" "EZH1cs" "EZH2" "H3K27me3" "IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103_heatmap_colorSmall3.pdf \
    --samplesLabel "SUZ12" "EZH1cs" "EZH2" "H3K27me3" "IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103_profile.pdf \
    --samplesLabel "SUZ12" "EZH1cs" "EZH2" "H3K27me3" "IGG" \
    --perGroup \
    --colors green purple blue red grey \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103.gz \
    -out output/deeptools/matrix_TSS_5kb_SUZ12_EZH1cs_EZH2_H3K27me3_IGG_raw_SUZ12Peaks_macs2q1.30103_heatmap_colorSmall4.pdf \
    --samplesLabel "SUZ12" "EZH1cs" "EZH2" "H3K27me3" "IGG" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10 \
    --zMax 10 10 10 80 80


