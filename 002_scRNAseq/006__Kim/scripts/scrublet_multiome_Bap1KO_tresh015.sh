#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00


python3 scripts/scrublet_doublets_histogram.py multiome_Bap1KO/outs/filtered_feature_bc_matrix output/doublets/multiome_Bap1KO_tresh015.tsv 0.15

