#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00





plotHeatmap -m output/deeptools/matrix_gene_1kb_DiffBind_TMM.gz \
    -out output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.png \
    --perGroup \
    --kmeans 6 \
    --samplesLabel "WT" "HET" "KO" \
    --outFileSortedRegions output/deeptools/matrix_gene_1kb_DiffBind_TMM_heatmap_kmeans6.txt



