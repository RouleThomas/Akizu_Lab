#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/gencode_upregulated_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-ESC_KO_vs_ESC_WT-H3K27me3.gtf meta/gencode_downregulated_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-ESC_KO_vs_ESC_WT-H3K27me3.gtf \
    -S output/bigwig_Ferguson/ESC_WT_H3K27me3_noXchr_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_KO_H3K27me3_noXchr_unique_norm99_initialBigwig_median.bw output/bigwig_Ferguson/ESC_OEKO_H3K27me3_noXchr_unique_norm99_initialBigwig_median.bw output/bigwig/ESC_WT_IGG_R1_noXchr.unique.dupmark.sorted.bw output/bigwig/ESC_WT_IGG_R2_noXchr.unique.dupmark.sorted.bw output/bigwig/ESC_WT_IGG_R3_noXchr.unique.dupmark.sorted.bw \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_heatmap.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_plotProfile1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colors black red blue grey grey grey \
    --perGroup \
    --plotWidth 8

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_plotProfile2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colors black red blue grey grey grey \
    --perGroup \
    --plotWidth 4

plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_plotProfile3.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colors black red blue grey grey grey \
    --perGroup \
    --plotWidth 4



# interactive


plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_heatmap1.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3




plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_heatmap2.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 2 2 2





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_heatmap3.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap Blues \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3







plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_heatmap4.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorList 'white,#e8f3ff,#cce6ff,#99ccff,#66b3ff,#3385cc,#0066cc,#004c99,#003366' \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_kmeans6_heatmap3.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap Blues \
    --kmeans 6 \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3



plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-gencode_q05fc058_promoter5-WTKOOEKO_H3K27me3_qval23merge100bpnoXchrthresh1-KO_vs_WT-WTKOOEKO-H3K27me3IGG_noSkip0_kmeans12_heatmap3.pdf \
    --samplesLabel "WT_H3K27me3" "KO_H3K27me3" "OEKO_H3K27me3" "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap Blues \
    --kmeans 12 \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3

