#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=8

computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../005__CutRun_NPC_PSC/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bw output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_threshold10.gz \
    -p 8 \
    --minThreshold 10




plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_threshold10.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_threshold10_heatmap.pdf \
    --samplesLabel "WT" "KO" "KOEF1aEZH1" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 
   # --zMin 0 0 0 --zMax 5 5 40 


plotProfile -m output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_threshold10.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_median_threshold10_profile.pdf \
    --samplesLabel "WT" "KO" "KOEF1aEZH1" \
    --refPointLabel "TSS" \
    --perGroup \
    --colors black red grey
  #  --yMin 0 0 0 --yMax 5 5 40 
    #  --colors lightblue purple pink \


