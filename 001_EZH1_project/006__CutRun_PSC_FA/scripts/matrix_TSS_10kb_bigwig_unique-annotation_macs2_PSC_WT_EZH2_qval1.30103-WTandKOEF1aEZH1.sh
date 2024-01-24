#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../006__CutRun_PSC_FA/output/macs2/broad_blacklist_qval1.30103/PSC_WT_EZH2_peaks.broadPeak \
    -S ../006__CutRun_PSC_FA/output/bigwig/PSC_WT_EZH2.unique.dupmark.sorted.bw ../006__CutRun_PSC_FA/output/bigwig/PSC_KOEF1aEZH1_EZH2.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_WT_EZH2_qval1.30103-WTandKOEF1aEZH1.gz \
    -p 8




plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_WT_EZH2_qval1.30103-WTandKOEF1aEZH1.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_WT_EZH2_qval1.30103-WTandKOEF1aEZH1_heatmap.pdf \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 
    # --zMin 0 0 0 --zMax 5 5 40 


plotProfile -m output/deeptools/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_WT_EZH2_qval1.30103-WTandKOEF1aEZH1.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_unique-annotation_macs2_PSC_WT_EZH2_qval1.30103-WTandKOEF1aEZH1_profile.pdf \
    --colors black darkblue \
    --samplesLabel "WT" "KOEF1aEZH1" \
    --refPointLabel "TSS" \
    --perGroup
   #  --yMin 0 0 0 --yMax 5 5 40 
