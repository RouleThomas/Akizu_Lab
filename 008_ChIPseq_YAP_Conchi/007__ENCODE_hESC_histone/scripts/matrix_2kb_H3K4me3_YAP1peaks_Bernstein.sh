#!/bin/bash
#SBATCH --mem=100G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 2000 -a 2000 \
    -R output/annotation_homer_hESC_WT_YAP1_R1_annot.bed \
    -S output/ENCODE/Bernstein_H3K4me1.bigWig output/ENCODE/Bernstein_H3K27ac.bigWig output/ENCODE/Bernstein_H3K36me3.bigWig output/ENCODE/Bernstein_H3K27me3.bigWig ../001__ChIPseq_V1/output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_2kb_H3K4me3_YAP1peaks_Bernstein.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_2kb_H3K4me3_YAP1peaks_Bernstein.gz \
    -out output/deeptools/matrix_2kb_H3K4me3_YAP1peaks_Bernstein_heatmap.pdf \
    --samplesLabel "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4


# interactive
plotHeatmap -m output/deeptools/matrix_2kb_H3K4me3_YAP1peaks_Bernstein.gz \
    -out output/deeptools/matrix_2kb_H3K4me3_YAP1peaks_Bernstein_heatmap1.pdf \
    --samplesLabel "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3

