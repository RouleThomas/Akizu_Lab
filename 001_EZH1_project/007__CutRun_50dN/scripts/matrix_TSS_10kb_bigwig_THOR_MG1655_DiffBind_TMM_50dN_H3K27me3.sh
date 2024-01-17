#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../005__CutRun_NPC_PSC/meta/ENCFF159KBI.gtf \
    -S output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1-rep0.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1-rep1.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2-rep0.bw output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2-rep1.bw output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2-rep0.bw output/THOR/THOR_50dN_H3K27me3_WTvsKOEF1aEZH1/50dNH3K27me3WTvsKOEF1aEZH1-s2-rep1.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3.gz \
    -p 6




plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_heatmap.pdf \
    --samplesLabel "WT_R1" "WT_R2" "KO_R1" "KO_R2" "KOEF1aEZH1_R1" "KOEF1aEZH1_R2" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 
   # --zMin 0 0 0 --zMax 5 5 40 


plotProfile -m output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_THOR_MG1655_DiffBind_TMM_50dN_H3K27me3_profile.pdf \
    --samplesLabel "WT_R1" "WT_R2" "KO_R1" "KO_R2" "KOEF1aEZH1_R1" "KOEF1aEZH1_R2" \
    --refPointLabel "TSS" \
    --perGroup \
    --colors black black red red grey grey
  #  --yMin 0 0 0 --yMax 5 5 40 
    #  --colors lightblue purple pink \