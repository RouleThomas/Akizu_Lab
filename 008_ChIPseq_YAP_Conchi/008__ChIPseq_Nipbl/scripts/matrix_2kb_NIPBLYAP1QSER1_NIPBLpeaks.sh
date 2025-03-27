#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 2000 -a 2000 \
    -R output/homer/Nipbl/peaks.bed \
    -S output/bigwig/Nipbl.unique.dupmark.sorted.bw ../003__ChIPseq_pluripotency/output/bigwig/hESC_WT_YAP1_R1.unique.dupmark.sorted.bw ../001__ChIPseq_V1/output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s1_median.bw  \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks.gz \
    -out output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks_heatmap.pdf \
    --samplesLabel "NIPBL" "YAP1" "QSER1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4


# interactive
plotHeatmap -m output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks.gz \
    -out output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks_heatmap1.pdf \
    --samplesLabel "NIPBL" "YAP1" "QSER1" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3



plotHeatmap -m output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks.gz \
    -out output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks_heatmap2.pdf \
    --samplesLabel "NIPBL" "YAP1" "QSER1" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3 \
    --zMax 10 4 10

