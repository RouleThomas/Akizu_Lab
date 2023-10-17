#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R meta/ENCFF159KBI_SUZ12-EZH2.gtf meta/ENCFF159KBI_SUZ12-noEZH2.gtf \
    -S output/THOR/THOR_NPC_SUZ12/NPCSUZ12-s1-rep0.bw output/THOR/THOR_NPC_EZH2/NPCEZH2-s1-rep0.bw output/THOR/THOR_NPC_H3K27me3/NPCH3K27me3-s1-rep0.bw output/THOR/THOR_NPC_H3K4me3/NPCH3K4me3-s1-rep0.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_THOR_SUZ12-EZH2_WT.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_THOR_SUZ12-EZH2_WT.gz \
    -out output/deeptools/matrix_TSS_10kb_THOR_SUZ12-EZH2_WT_heatmap.pdf \
    --samplesLabel "SUZ12" "EZH2" "H3K27me3" "H3K4me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 \
    --zMin 0 0 0 0 --zMax 20 20 100 100