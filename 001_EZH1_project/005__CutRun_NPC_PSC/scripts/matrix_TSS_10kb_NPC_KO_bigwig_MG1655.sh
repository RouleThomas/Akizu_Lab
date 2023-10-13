#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/bigwig_MG1655_DiffBind_TMM/NPC_KO_EZH2.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_SUZ12.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27me3.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3.unique.dupmark.sorted.bw output/bigwig_MG1655_DiffBind_TMM/NPC_KO_IGG.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_NPC_KO_bigwig_MG1655.gz \
    -p 8


plotHeatmap -m output/deeptools/matrix_TSS_10kb_NPC_KO_bigwig_MG1655.gz \
    -out output/deeptools/matrix_TSS_10kb_NPC_KO_bigwig_MG1655_heatmap.png \
    --samplesLabel "EZH2" "SUZ12" "H3K27me3" "H3K4me3" "IGG" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 4



plotProfile -m output/deeptools/matrix_TSS_10kb_NPC_KO_bigwig_MG1655.gz \
    -out output/deeptools/matrix_TSS_10kb_NPC_KO_bigwig_MG1655_profile.pdf \
    --samplesLabel "EZH2" "SUZ12" "H3K27me3" "H3K4me3" "IGG" \
    --perGroup \
    --colors blue green purple orange gray \
    --refPointLabel "TSS" \
    -T "H3K27me3 read density" \
    -z ""


