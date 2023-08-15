#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix scale-regions \
    -b 5000 -a 5000 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_THOR_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOR_allGenes.gz \
    -out output/deeptools/matrix_TSS_5kb_THOR_allGenes_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList white,grey,black \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15

