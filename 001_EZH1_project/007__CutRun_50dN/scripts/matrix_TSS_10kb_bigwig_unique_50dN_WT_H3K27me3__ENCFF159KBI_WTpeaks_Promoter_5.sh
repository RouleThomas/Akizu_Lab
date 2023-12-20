#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../003__CutRun/meta/ENCFF159KBI_WTpeaks_Promoter_5.gtf \
    -S output/bigwig/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bw output/bigwig/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bw output/bigwig/50dN_WTQ731E_H3K27me3_R3.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_bigwig_unique_50dN_WT_H3K27me3__ENCFF159KBI_WTpeaks_Promoter_5.gz \
    -p 6




plotHeatmap -m output/deeptools/matrix_TSS_10kb_bigwig_unique_50dN_WT_H3K27me3__ENCFF159KBI_WTpeaks_Promoter_5.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_unique_50dN_WT_H3K27me3__ENCFF159KBI_WTpeaks_Promoter_5_heatmap2.pdf \
    --samplesLabel "H3K27me3_R1" "H3K27me3_R2" "H3K27me3_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 15 \
    --heatmapWidth 2 
   # --zMin 0 0 0 --zMax 5 5 40 


plotProfile -m output/deeptools/matrix_TSS_10kb_bigwig_unique_50dN_WT_H3K27me3__ENCFF159KBI_WTpeaks_Promoter_5.gz \
    -out output/deeptools/matrix_TSS_10kb_bigwig_unique_50dN_WT_H3K27me3__ENCFF159KBI_WTpeaks_Promoter_5_profile.pdf \
    --samplesLabel "H3K27me3_R1" "H3K27me3_R2" "H3K27me3_R3" \
    --refPointLabel "TSS" \
    --perGroup
  #  --yMin 0 0 0 --yMax 5 5 40 
    #  --colors lightblue purple pink \