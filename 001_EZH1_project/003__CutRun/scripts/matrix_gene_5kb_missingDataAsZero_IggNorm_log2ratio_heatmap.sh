#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00


plotHeatmap -m output/deeptools/matrix_gene_5kb_missingDataAsZero_IggNorm_log2ratio.gz \
    -out output/deeptools/matrix_gene_5kb_missingDataAsZero_IggNorm_log2ratio_heatmap.png \
    --colorMap RdBu \
    --whatToShow 'heatmap and colorbar'



