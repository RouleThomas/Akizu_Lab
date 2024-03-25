#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



computeMatrix reference-point --referencePoint center \
    -b 10000 -a 10000 \
    -R ../../001_EZH1_Project/007__CutRun_50dN/output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20_positive.bed ../../001_EZH1_Project/007__CutRun_50dN/output/THOR/THOR_50dN_H3K27me3_WTvsKO/THOR_qval20_negative.bed \
    -S ../../001_EZH1_Project/007__CutRun_50dN/output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s1_median.bw ../../001_EZH1_Project/007__CutRun_50dN/output/THOR/THOR_50dN_H3K27me3_WTvsKO/50dNH3K27me3WTvsKO-s2_median.bw output/THOR/THOR_EZH2inh_H3K27me3/EZH2inhH3K27me3-s1_median.bw output/THOR/THOR_EZH2inh_H3K27me3/EZH2inhH3K27me3-s2_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks_heatmap.pdf \
    --samplesLabel "WT" "KO" "DMSO" "EZH2inh" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2


plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks_heatmap1.pdf \
    --samplesLabel "WT" "KO" "DMSO" "EZH2inh" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue \
    --colorNumber 2 \
    --zMax 120 120 10 10 10 10




plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks_heatmap2.pdf \
    --samplesLabel "WT" "KO" "DMSO" "EZH2inh" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList grey,blue,blue,blue \
    --colorNumber 4 \
    --zMax 120 120 10 10 10 10



plotHeatmap -m output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks_heatmap3.pdf \
    --samplesLabel "WT" "KO" "DMSO" "EZH2inh" \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --colorList White,White,LightGrey,"#c5c9c7",DarkGray,Gray,DimGray,"#3c4142","#1b2431",Black \
    --colorNumber 10 \
    --zMax 120 120 10 10 10 10



plotProfile -m output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks.gz \
    -out output/deeptools/matrix_TSS_10kb_H3K27me3_CutRun007_CiceriEpiInh_007H3K72me3GainLost_007THOR_q20_peaks_profile.pdf \
    --samplesLabel "WT" "KO" "DMSO" "EZH2inh" \
    --colors black red darkgrey darkred \
    -T "Read density" \
    --numPlotsPerRow 2



