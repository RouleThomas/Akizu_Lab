#!/bin/bash
#SBATCH --mem=350G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7


computeMatrix scale-regions \
    -b 1000 -a 300 \
    -R ../../Master/meta/ENCFF159KBI.gtf \
    -S ../001__RNAseq/output/bigwig_hg38/8wN_WT_median.bw output/THOR/THOR_WTvsHET_unique_Keepdup/WTvsHETuniqueKeepdup-s1_median.bw \
    --skipZeros \
    --missingDataAsZero \
    --blackListFileName ../../Master/meta/hg38-blacklist.v2.bed \
    -o output/deeptools/matrix_gene_1kb300bp_TPM_THOR_WT_allGenes.gz \
    -p 6



plotHeatmap -m output/deeptools/matrix_gene_1kb300bp_TPM_THOR_WT_allGenes.gz \
    -out output/deeptools/matrix_gene_1kb300bp_TPM_THOR_WT_allGenes_heatmap.pdf \
    --samplesLabel "RNA" "H3K27me3" \
    --colorMap bwr \
    --whatToShow 'heatmap and colorbar' \
    --heatmapHeight 10 \
    --heatmapWidth 2 \
    --yMax 0.1 10

