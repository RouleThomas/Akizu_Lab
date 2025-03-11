#!/bin/bash
#SBATCH --mem=200G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6




computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI__NPC_WTKO_H3K27me3_SICER2window200gap600fdr05evalue50000_gain_annot_promoterAnd5.gtf meta/ENCFF159KBI__NPC_WTKO_H3K27me3_SICER2window200gap600fdr05evalue50000_lost_annot_promoterAnd5.gtf \
    -S output/bigwig_Ferguson/NPC_WT_H3K27me3_unique_norm99_median.bw output/bigwig_Ferguson/NPC_KO_H3K27me3_unique_norm99_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-gene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-gene_heatmap.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 10 10




plotHeatmap -m output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-gene_heatmap1.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 5 5



plotHeatmap -m output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-gene.gz \
    -out output/deeptools/matrix_TSS_5kb-SICER2-SICER2window200gap600fdr05evalue50000-gene_heatmap2.pdf \
    --samplesLabel "WT_LocalMaxima" "KO_LocalMaxima" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 8 \
    --heatmapWidth 2 \
    --zMax 2 2