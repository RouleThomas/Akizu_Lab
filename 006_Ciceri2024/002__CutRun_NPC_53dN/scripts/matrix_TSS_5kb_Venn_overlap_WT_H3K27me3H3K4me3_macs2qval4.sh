#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K4me3onlyqval4.gtf meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3andH3K4me3qval4.gtf meta/ENCFF159KBI_Venn_overlap_WT_H3K27me3H3K4me3__H3K27me3onlyqval4.gtf \
    -S output/bigwig/53dN_WT_H3K4me3_median.bw output/bigwig/53dN_WT_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4.gz \
    -p 6




plotHeatmap -m output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4_heatmap_colorSmall.pdf \
    --samplesLabel "H3K4me3" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4_heatmap_colorSmall2.pdf \
    --samplesLabel "H3K4me3" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 50 50

plotHeatmap -m output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4_heatmap_colorSmall3.pdf \
    --samplesLabel "H3K4me3" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --zMax 30 30


plotProfile -m output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4.gz \
    -out output/deeptools/matrix_TSS_5kb_Venn_overlap_WT_H3K27me3H3K4me3_macs2qval4_profile.pdf \
    --samplesLabel "H3K4me3" "H3K27me3" \
    --perGroup \
    --colors orange purple \
    --refPointLabel "TSS"



