#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_WTpeaks_Promoter_5.gtf  \
    -S output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s2_median.bw output/THOR/THOR_WTvsKO_unique_Keepdup/WTvsKOuniqueKeepdup-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_THOR_WTpeaksGene.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOR_WTpeaksGene.gz \
    -out output/deeptools/matrix_TSS_5kb_THOR_WTpeaksGene_heatmap.png \
    --samplesLabel "WT" "HET" "KO" \
    --colorList white,grey,black \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --refPointLabel 0
