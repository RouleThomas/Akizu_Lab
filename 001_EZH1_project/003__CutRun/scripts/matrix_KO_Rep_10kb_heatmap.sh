#!/bin/bash
#SBATCH --mem=400G
#SBATCH --time=200:00:00



plotHeatmap -m output/deeptools/matrix_KO_Rep_10kb.gz \
-out output/deeptools/matrix_KO_Rep_10kb_heatmap.png \
--colorMap RdBu \
--whatToShow 'heatmap and colorbar'



