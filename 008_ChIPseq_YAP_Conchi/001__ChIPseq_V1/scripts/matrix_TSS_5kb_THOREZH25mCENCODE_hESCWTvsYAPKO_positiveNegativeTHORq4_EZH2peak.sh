#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_positive.bed output/THOR/THOR_hESC_EZH2_WTvsYAPKO/THOR_qval4_negative.bed \
    -S output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s2_median.bw ../004__WGBS_Dixon2021/output/ENCODE/ENCFF969PDB.bigWig ../004__WGBS_Dixon2021/output/ENCODE/ENCFF423PGW.bigWig ../004__WGBS_Dixon2021/output/ENCODE/ENCFF725YJG.bigWig ../004__WGBS_Dixon2021/output/ENCODE/ENCFF040LKO.bigWig \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak_heatmap_colorSmall1.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Ecker_R1" "Ecker_R2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 40 40 0.1 0.1 1 1


plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak_heatmap_colorSmall.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Ecker_R1" "Ecker_R2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak_heatmap_colorSmall2.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Ecker_R1" "Ecker_R2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak_heatmap_colorSmall3.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Ecker_R1" "Ecker_R2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak_profile.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Ecker_R1" "Ecker_R2" "Myers_R1" "Myers_R2" \
    --perGroup \
    --colors black darkred grey grey grey grey \
    --refPointLabel "center" \
    -T "Read density" \
    --plotWidth 10


plotHeatmap -m output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak.gz \
    -out output/deeptools/matrix_TSS_5kb_THOREZH25mCENCODE_hESCWTvsYAPKO_positiveNegativeTHORq4_EZH2peak_heatmap_colorSmall4.pdf \
    --samplesLabel "WT_EZH2" "YAPKO_EZH2" "Ecker_R1" "Ecker_R2" "Myers_R1" "Myers_R2" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



