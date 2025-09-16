#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix reference-point --referencePoint center \
    -b 5000 -a 5000 \
    -R output/macs2/broad/broad_blacklist_qval2.30103/ESC_WTKOOEKO_H3K27me3_noXchr_pool_peaks.sorted.merge100bp.bed \
    -S output/bigwig/ESC_WT_IGG_R1.unique.dupmark.sorted.bw output/bigwig/ESC_WT_IGG_R2.unique.dupmark.sorted.bw output/bigwig/ESC_WT_IGG_R3.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3_heatmap.pdf \
    --samplesLabel "WT_IGG_R1" "WT_IGG_R2" "WT_IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3_plotProfile1.pdf \
    --samplesLabel "WT_IGG_R1" "WT_IGG_R2" "WT_IGG_R3" \
    --colors grey grey grey \
    --perGroup \
    --plotWidth 8



# interactive


plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3_heatmap1.pdf \
    --samplesLabel "WT_IGG_R1" "WT_IGG_R2" "WT_IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 3 3 3




plotHeatmap -m output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3.gz \
    -out output/deeptools/matrix_PEAK_5kb-macs2broad_WTKOOEKOconsensus_H3K27me3poolqval23merge100bp-WTKOOEKO-IGGR1R2R3_heatmap2.pdf \
    --samplesLabel "WT_IGG_R1" "WT_IGG_R2" "WT_IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 1 1 1 


