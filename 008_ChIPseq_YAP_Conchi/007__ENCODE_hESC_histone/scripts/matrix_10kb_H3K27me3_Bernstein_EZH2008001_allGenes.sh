#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI.gtf \
    -S ../001__ChIPseq_V1/output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw output/ENCODE/Bernstein_H3K27me3.bigWig \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_10kb_H3K27me3_Bernstein_EZH2008001_allGenes.gz \
    -p 6




plotHeatmap -m output/deeptools/matrix_10kb_H3K27me3_Bernstein_EZH2008001_allGenes.gz \
    -out output/deeptools/matrix_10kb_H3K27me3_Bernstein_EZH2008001_allGenes_heatmap.pdf \
    --samplesLabel "EZH2" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4



# interactive
plotHeatmap -m output/deeptools/matrix_10kb_H3K27me3_Bernstein_EZH2008001_allGenes.gz \
    -out output/deeptools/matrix_10kb_H3K27me3_Bernstein_EZH2008001_allGenes_heatmap1.pdf \
    --samplesLabel "EZH2" "H3K27me3" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 5


