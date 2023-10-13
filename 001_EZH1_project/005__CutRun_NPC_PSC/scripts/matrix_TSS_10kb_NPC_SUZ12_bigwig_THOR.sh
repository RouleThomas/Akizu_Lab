#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint TSS \
    -b 10000 -a 10000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s1-rep0.bw output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_NPC_SUZ12_bigwig_THOR.gz \
    -p 8


plotHeatmap -m output/deeptools/matrix_TSS_10kb_NPC_SUZ12_bigwig_THOR.gz \
    -out output/deeptools/matrix_TSS_10kb_NPC_SUZ12_bigwig_THOR_heatmap.png \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 4



plotProfile -m output/deeptools/matrix_TSS_10kb_NPC_SUZ12_bigwig_THOR.gz \
    -out output/deeptools/matrix_TSS_10kb_NPC_SUZ12_bigwig_THOR_profile.pdf \
    --samplesLabel "WT" "KO" \
    --perGroup \
    --colors black red \
    --refPointLabel "TSS" \
    -T "SUZ12 read density" \
    -z ""


