#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../005__CutRun_NPC_PSC/meta/ENCFF159KBI.gtf \
    -S output/bigwig_MG1655_DiffBind_TMM/50dN_WTQ731E_H3K27me3_median.bw output/bigwig_MG1655_DiffBind_TMM/50dN_KO_H3K27me3_median.bw output/bigwig_MG1655_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_H3K27me3_median.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_H3K27me3_median.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_H3K27me3_median_heatmap.pdf \
    --samplesLabel "WT" "KO" "KOEF1aEZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 
   # --zMin 0 0 0 --zMax 5 5 40 


plotProfile -m output/deeptools/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_H3K27me3_median.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_unique_MG1655_DiffBind_TMM_50dN_H3K27me3_median_profile.pdf \
    --samplesLabel "WT" "KO" "KOEF1aEZH1" \
    --refPointLabel "TSS" \
    --perGroup \
    --colors black red grey
  #  --yMin 0 0 0 --yMax 5 5 40 
    #   \