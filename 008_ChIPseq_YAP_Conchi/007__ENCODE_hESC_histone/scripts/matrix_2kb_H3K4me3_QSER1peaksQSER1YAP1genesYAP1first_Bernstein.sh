#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 2000 -a 2000 \
    -R ../001__ChIPseq_V1/output/homer/QSER1_YAP1_1192genes_ChIPseeker_noHeader.bed \
    -S ../003__ChIPseq_pluripotency/output/bigwig/hESC_WT_YAP1_R1.unique.dupmark.sorted.bw ../001__ChIPseq_V1/output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s1_median.bw output/ENCODE/Bernstein_H3K4me1.bigWig output/ENCODE/Bernstein_H3K27ac.bigWig output/ENCODE/Bernstein_H3K36me3.bigWig output/ENCODE/Bernstein_H3K27me3.bigWig ../001__ChIPseq_V1/output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw ../008__ChIPseq_Nipbl/output/bigwig/Nipbl.unique.dupmark.sorted.bw ../008__ChIPseq_Nipbl/output/bigwig/input.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein.gz \
    -out output/deeptools/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein_heatmap.pdf \
    --samplesLabel "YAP1" "QSER1" "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" "NIPBL" "input" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4


# interactive
plotHeatmap -m output/deeptools/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein.gz \
    -out output/deeptools/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein_heatmap1.pdf \
    --samplesLabel "YAP1" "QSER1" "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" "NIPBL" "input" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3

plotHeatmap -m output/deeptools/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein.gz \
    -out output/deeptools/matrix_2kb_H3K4me3_QSER1peaksQSER1YAP1genesYAP1first_Bernstein_heatmap2.pdf \
    --samplesLabel "YAP1" "QSER1" "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" "NIPBL" "input" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3 \
    --zMax 5 20 15 15 15 15 15 15 15 15 
