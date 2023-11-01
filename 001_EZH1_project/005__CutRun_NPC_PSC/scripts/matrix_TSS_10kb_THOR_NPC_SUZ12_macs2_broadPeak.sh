#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/macs2/broad_blacklist_qval1.30103/NPC_WT_SUZ12_peaks.broadPeak \
    -S output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s1-rep0.bw output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_THOR_NPC_SUZ12_macs2_broadPeak.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_NPC_SUZ12_macs2_broadPeak.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_NPC_SUZ12_macs2_broadPeak_heatmap.pdf \
    --samplesLabel "WT" "KO" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2

