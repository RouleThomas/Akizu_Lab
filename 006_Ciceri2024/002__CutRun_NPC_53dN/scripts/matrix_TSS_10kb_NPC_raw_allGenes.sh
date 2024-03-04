#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig/NPC_WT_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27ac_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27ac_R2.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K4me3_R2.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K9me3_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K9me3_R2.unique.dupmark.sorted.bw output/bigwig/NPC_WT_IGG_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_IGG_R2.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_NPC_raw_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_NPC_raw_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_NPC_raw_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "H3K27me3_R1" "H3K27me3_R2" "H3K27ac_R1" "H3K27ac_R2" "H3K4me3_R1" "H3K4me3_R2" "H3K9me3_R1" "H3K9me3_R2" "IGG_R1" "IGG_R2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2

plotProfile -m output/deeptools/matrix_TSS_10kb_NPC_raw_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_NPC_raw_allGenes_profile.pdf \
    --samplesLabel "H3K27me3_R1" "H3K27me3_R2" "H3K27ac_R1" "H3K27ac_R2" "H3K4me3_R1" "H3K4me3_R2" "H3K9me3_R1" "H3K9me3_R2" "IGG_R1" "IGG_R2" \
    --perGroup \
    --colors green green red red blue blue purple purple grey grey \
    --refPointLabel "TSS" \
    -T "Read density" \
    -z ""



