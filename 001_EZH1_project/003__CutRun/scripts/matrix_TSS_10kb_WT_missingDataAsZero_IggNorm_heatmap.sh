#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00




plotHeatmap -m output/deeptools/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm.gz \
    -out output/deeptools/matrix_TSS_10kb_WT_missingDataAsZero_IggNorm_heatmap.png \
    --colorMap RdBu \
    --whatToShow 'heatmap and colorbar'


