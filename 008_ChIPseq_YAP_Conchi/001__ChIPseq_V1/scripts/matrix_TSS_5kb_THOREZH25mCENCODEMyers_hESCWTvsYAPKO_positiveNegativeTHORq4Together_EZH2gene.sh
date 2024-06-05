#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_THORq4_EZH2_posneg_annot_promoterAnd5.gtf \
    -S output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2_median.bw ../004__WGBS_Dixon2021/output/ENCODE/ENCFF725YJG.bigWig ../004__WGBS_Dixon2021/output/ENCODE/ENCFF040LKO.bigWig \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene_heatmap_colorSmall1.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 40 40 1 1


plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene_heatmap_colorSmall.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene_heatmap_colorSmall3.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene_profile.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Myers_R1" "Myers_R2" \
    --colors black grey \
    --refPointLabel "center" \
    -T "Read density" \
    --plotWidth 10 \
    --yMax 1


plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODEMyers_hESCWTvsYAPKO_positiveNegativeTHORq4Together_EZH2gene_heatmap_colorSmall4.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



