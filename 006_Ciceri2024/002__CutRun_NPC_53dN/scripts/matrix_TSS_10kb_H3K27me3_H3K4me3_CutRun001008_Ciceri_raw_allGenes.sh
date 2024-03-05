#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K27me3_005.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K27me3_008.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K27me3_R2.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K4me3_005.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_H3K4me3_008.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K4me3_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_H3K4me3_R2.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_IGG_005.unique.dupmark.sorted.bw ../../001_EZH1_Project/009__integration_NPC_WTKO_K27me3K4me3_005_008/output/bigwig/NPC_WT_IGG_008.unique.dupmark.sorted.bw output/bigwig/NPC_WT_IGG_R1.unique.dupmark.sorted.bw output/bigwig/NPC_WT_IGG_R2.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_H3K4me3_CutRun001008_Ciceri_raw_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_H3K4me3_CutRun001008_Ciceri_raw_allGenes.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_H3K4me3_CutRun001008_Ciceri_raw_allGenes_heatmap_colorSmall.pdf \
    --samplesLabel "Akizu_H3K27me3_R1" "Akizu_H3K27me3_R2" "Ciceri_H3K27me3_R1" "Ciceri_H3K27me3_R2" "Akizu_H3K4me3_R1" "Akizu_H3K4me3_R2" "Ciceri_H3K4me3_R1" "Ciceri_H3K4me3_R2" "Akizu_IGG_R1" "Akizu_IGG_R2" "Ciceri_IGG_R1" "Ciceri_IGG_R2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2




