#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WT_SUZ12_macs2_broadPeak_heatmap_kmeans25_cluster25.bed \
    -S output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s1-rep0.bw output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s2-rep0.bw output/THOR/THOR_NPC_EZH2/NPCEZH2-s1-rep0.bw output/THOR/THOR_NPC_EZH2/NPCEZH2-s2-rep0.bw output/THOR/THOR_NPC_H3K27me3/NPCH3K27me3-s1-rep0.bw output/THOR/THOR_NPC_H3K27me3/NPCH3K27me3-s2-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WTandKO_SUZ12_macs2_broadPeak_cl25.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WTandKO_SUZ12_macs2_broadPeak_cl25.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12_EZH2_WTandKO_SUZ12_macs2_broadPeak_cl25_heatmap.pdf \
    --samplesLabel "SUZ12_WT" "SUZ12_KO" "EZH2_WT" "EZH2_KO" "H3K27me3_WT" "H3K27me3_KO" \
    --colorMap bwr \
    --whatToShow 'plot, heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --zMin 0 0 0 0 0 0 --zMax 30 30 30 30 50 50

