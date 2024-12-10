#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_homer_hESC_WT_EZH2_pool.gtf \
    -S ../001__ChIPseq_V1/output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw output/ENCODE/Bernstein_H3K27me3.bigWig \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_5kb_H3K27me3_Bernstein_EZH2008001_EZH2target.gz \
    -p 6




plotHeatmap -m output/deeptools/matrix_5kb_H3K27me3_Bernstein_EZH2008001_EZH2target.gz \
    -out output/deeptools/matrix_5kb_H3K27me3_Bernstein_EZH2008001_EZH2target_heatmap.pdf \
    --samplesLabel "EZH2" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 4



# interactive
plotHeatmap -m output/deeptools/matrix_5kb_H3K27me3_Bernstein_EZH2008001_EZH2target.gz \
    -out output/deeptools/matrix_5kb_H3K27me3_Bernstein_EZH2008001_EZH2target_heatmap1.pdf \
    --samplesLabel "EZH2" "H3K27me3" \
    --colorList 'black, yellow' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 5


