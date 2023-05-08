#!/bin/bash
#SBATCH --mem=300G
#SBATCH --time=200:00:00



plotHeatmap -m output/deeptools/matrix_TSS_10kb_WT_missingDataAsZero.gz \
-out output/deeptools/matrix_TSS_10kb_WT_missingDataAsZero_heatmap.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar'


