#!/bin/bash
#SBATCH --mem=399G
#SBATCH --time=100:00:00


python3 scripts/scrublet_doublets.py WT_p180_CB_Rep1/outs/filtered_feature_bc_matrix output/doublets/WT_p180_CB_Rep1.tsv
python3 scripts/scrublet_doublets.py WT_p180_CB_Rep2/outs/filtered_feature_bc_matrix output/doublets/WT_p180_CB_Rep2.tsv
python3 scripts/scrublet_doublets.py WT_p180_CB_Rep3/outs/filtered_feature_bc_matrix output/doublets/WT_p180_CB_Rep3.tsv
python3 scripts/scrublet_doublets.py WT_p180_CX_Rep1/outs/filtered_feature_bc_matrix output/doublets/WT_p180_CX_Rep1.tsv
python3 scripts/scrublet_doublets.py WT_p180_CX_Rep2/outs/filtered_feature_bc_matrix output/doublets/WT_p180_CX_Rep2.tsv
python3 scripts/scrublet_doublets.py WT_p180_CX_Rep3/outs/filtered_feature_bc_matrix output/doublets/WT_p180_CX_Rep3.tsv


