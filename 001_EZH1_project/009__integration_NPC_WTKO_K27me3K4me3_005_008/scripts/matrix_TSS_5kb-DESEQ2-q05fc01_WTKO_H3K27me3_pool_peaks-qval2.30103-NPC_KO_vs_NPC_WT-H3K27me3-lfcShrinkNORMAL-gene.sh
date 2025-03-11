#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI__NPC_WTKO_H3K27me3_upregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.gtf meta/ENCFF159KBI__NPC_WTKO_H3K27me3_downregulated_q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL.gtf \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb-DESEQ2-q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL-gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb-DESEQ2-q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-DESEQ2-q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL-gene_heatmap.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 10 10



plotHeatmap -m output/deeptools/matrix_TSS_5kb-DESEQ2-q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-DESEQ2-q05fc01_WTKO_H3K27me3_pool_peaks-qval2.30103-NPC_KO_vs_NPC_WT-H3K27me3-lfcShrinkNORMAL-gene_heatmap1.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 2 2


