#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6


computeMatrix reference-point --referencePoint center \
    -b 2000 -a 2000 \
    -R output/annotation_homer_hESC_WT_QSER1_pool_annot.bed \
    -S ../001__ChIPseq_V1/output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s1_median.bw output/ENCODE/Ren_H3K4me1.bigWig output/ENCODE/Ren_H3K27ac.bigWig output/ENCODE/Ren_H3K36me3.bigWig output/ENCODE/Ren_H3K27me3.bigWig ../001__ChIPseq_V1/output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw ../008__ChIPseq_Nipbl/output/bigwig/Nipbl.unique.dupmark.sorted.bw ../008__ChIPseq_Nipbl/output/bigwig/input.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_2kb_QSER1peaks_Ren_2input.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_2kb_QSER1peaks_Ren_2input.gz \
    -out output/deeptools/matrix_2kb_QSER1peaks_Ren_2input_heatmap.pdf \
    --samplesLabel "QSER1" "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" "NIPBL" "input" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4


# interactive
plotHeatmap -m output/deeptools/matrix_2kb_QSER1peaks_Ren_2input.gz \
    -out output/deeptools/matrix_2kb_QSER1peaks_Ren_2input_heatmap1.pdf \
    --samplesLabel "QSER1" "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" "NIPBL" "input" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3

plotHeatmap -m output/deeptools/matrix_2kb_QSER1peaks_Ren_2input.gz \
    -out output/deeptools/matrix_2kb_QSER1peaks_Ren_2input_heatmap2.pdf \
    --samplesLabel "QSER1" "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" "NIPBL" "input" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3 \
    --zMax 25 10 10 5 20 15 20 20


plotHeatmap -m output/deeptools/matrix_2kb_QSER1peaks_Ren_2input.gz \
    -out output/deeptools/matrix_2kb_QSER1peaks_Ren_2input_heatmap3.pdf \
    --samplesLabel "QSER1" "H3K4me1" "H3K27ac" "H3K36me3" "H3K27me3" "EZH2" "NIPBL" "input" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 3 \
    --zMax 25 8 8 4 15 10 10 10

