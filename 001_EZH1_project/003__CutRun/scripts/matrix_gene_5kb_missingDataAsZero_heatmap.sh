#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00




plotHeatmap -m output/deeptools/matrix_gene_5kb_missingDataAsZero.gz \
    -out output/deeptools/matrix_gene_5kb_missingDataAsZero_heatmap.png \
    --colorMap RdBu \
    --whatToShow 'heatmap and colorbar'



