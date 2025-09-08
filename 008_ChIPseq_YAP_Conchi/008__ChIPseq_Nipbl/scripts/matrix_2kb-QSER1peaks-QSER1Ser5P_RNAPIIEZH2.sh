#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 2000 -a 2000 \
    -R ../007__ENCODE_hESC_histone/output/annotation_homer_hESC_WT_QSER1_pool_annot.bed \
    -S ../001__ChIPseq_V1/output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s1_median.bw output/bigwig/Ser5P_RNAPII.unique.dupmark.sorted.bw ../001__ChIPseq_V1/output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2.gz \
    -out output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2_heatmap.pdf \
    --samplesLabel "QSER1" "Ser5P_RNAPII" "EZH2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4


# interactive
plotHeatmap -m output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2.gz \
    -out output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2_heatmap1.pdf \
    --samplesLabel "QSER1" "Ser5P_RNAPII" "EZH2" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3



plotHeatmap -m output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2.gz \
    -out output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2_heatmap2.pdf \
    --samplesLabel "QSER1" "Ser5P_RNAPII" "EZH2" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3 \
    --zMax 10 5 10


plotHeatmap -m output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2.gz \
    -out output/deeptools/matrix_2kb-QSER1peaks-QSER1Ser5P_RNAPIIEZH2_heatmap3.pdf \
    --samplesLabel "QSER1" "Ser5P_RNAPII" "EZH2" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3 \
    --zMax 10 5 5