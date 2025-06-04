#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 2000 -a 2000 \
    -R output/homer/Nipbl/peaks.bed \
    -S ../001__ChIPseq_V1/output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s1_median.bw output/bigwig/Nipbl.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL.gz \
    -out output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL_heatmap.pdf \
    --samplesLabel "QSER1" "NIPBL" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


# interactive
plotHeatmap -m output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL.gz \
    -out output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL_heatmap1.pdf \
    --samplesLabel "QSER1" "NIPBL" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotHeatmap -m output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL.gz \
    -out output/deeptools/matrix_2kb_NIPBLYAP1QSER1_NIPBLpeaks-QSER1NIPBL_heatmap2.pdf \
    --samplesLabel "QSER1" "NIPBL" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 10 10

