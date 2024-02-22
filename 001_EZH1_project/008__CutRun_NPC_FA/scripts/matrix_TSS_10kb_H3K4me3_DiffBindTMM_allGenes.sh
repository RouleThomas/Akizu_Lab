#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K4me3.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K4me3_DiffBindTMM_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K4me3_DiffBindTMM_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K4me3_DiffBindTMM_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

plotProfile -m output/deeptools/matrix_TSS_10kb_H3K4me3_DiffBindTMM_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K4me3_DiffBindTMM_allGenes_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "H3K4me3 read density" \
    -z ""



