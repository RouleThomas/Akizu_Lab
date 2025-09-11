#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=100:00:00
#SBATCH --cpus-per-task=6

computeMatrix scale-regions \
    -b 250 -a 100 \
    -R meta/ENCFF159KBI_IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO.gtf \
    -S output/bigwig/ESC_WT_IGG_R1.unique.dupmark.sorted.bw output/bigwig/ESC_WT_IGG_R2.unique.dupmark.sorted.bw output/bigwig/ESC_WT_IGG_R3.unique.dupmark.sorted.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3.gz \
    -p 6





plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3_heatmap.pdf \
    --samplesLabel "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2



plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3_plotProfile1.pdf \
    --samplesLabel "IGG_R1" "IGG_R2" "IGG_R3" \
    --colors black grey lightgrey \
    --perGroup \
    --plotWidth 4


plotProfile -m output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3_plotProfile2.pdf \
    --samplesLabel "IGG_R1" "IGG_R2" "IGG_R3" \
    --colors black grey lightgrey \
    --perGroup \
    --plotWidth 8


# interactive


plotHeatmap -m output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3.gz \
    -out output/deeptools/matrix_GENETSSTES_250bp100bp-IsoformSwitchAnalyzeR_kallisto_significant_isoforms_dIF01qval05__KO-WTKOOEKO-rawIGGR1R2R3_heatmap1.pdf \
    --samplesLabel "IGG_R1" "IGG_R2" "IGG_R3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 7 \
    --heatmapWidth 2 \
    --zMax 1 1 1




