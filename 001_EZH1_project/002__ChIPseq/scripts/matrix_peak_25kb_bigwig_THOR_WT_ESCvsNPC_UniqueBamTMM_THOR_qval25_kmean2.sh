#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


plotHeatmap -m output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25.gz \
    -out output/deeptools/matrix_peak_25kb_bigwig_THOR_WT_ESCvsNPC_UniqueBamTMM_THOR_qval25_kmean2_heatmap.png \
    --colorMap RdBu \
    --refPointLabel center \
    --kmeans 2


