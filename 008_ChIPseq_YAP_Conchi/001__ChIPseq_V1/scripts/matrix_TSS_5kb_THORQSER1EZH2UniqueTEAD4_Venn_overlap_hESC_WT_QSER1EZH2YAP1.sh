#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6          # Number of CPU cores per task



computeMatrix reference-point --referencePoint TSS \
    -b 5000 -a 5000 \
    -R meta/ENCFF159KBI_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gtf \
    -S output/THOR/THOR_hESC_QSER1_WTvsYAPKO/hESCQSER1WTvsYAPKO-s1_median.bw output/THOR/THOR_hESC_EZH2_WTvsYAPKO/hESCEZH2WTvsYAPKO-s1_median.bw ../003__ChIPseq_pluripotency/output/bigwig/hESC_WT_YAP1_R1.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
    -p 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
    -out output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1_heatmap_colorSmall1.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorMap bwr \
    --zMax 16 15 4


plotHeatmap -m output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
    -out output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1_heatmap_colorSmall.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
    -out output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1_heatmap_colorSmall2.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,darkblue,darkblue,darkblue \
    --colorNumber 6

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
    -out output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1_heatmap_colorSmall3.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4


plotProfile -m output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
    -out output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1_profile.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" \
    --perGroup \
    --colors blue red yellow \
    --refPointLabel "TSS" \
    -T "Read density" \
    -z ""

plotHeatmap -m output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
    -out output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1_heatmap_colorSmall4.pdf \
    --samplesLabel "QSER1" "EZH2" "YAP1" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10



plotProfile -m output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1.gz \
      --perGroup \
      --plotType heatmap \
      -out output/deeptools/matrix_TSS_5kb_THORQSER1EZH2UniqueTEAD4_Venn_overlap_hESC_WT_QSER1EZH2YAP1_profileHeatmap.pdf \
      --plotHeight 5 \
      --plotWidth 15



