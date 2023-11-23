#!/bin/bash
#SBATCH --mem=20G
#SBATCH --time=72:00:00


x=("S_CB_KO1" "S_CB_KO2" "S_CB_KO3" "S_CB_WT1" "S_CB_WT2" "S_CB_WT3" "S_CX_KO1" "S_CX_KO2" "S_CX_KO3" "S_CX_WT1" "S_CX_WT2" "S_CX_WT3" "171HetCB" "174MTCB" "175HetCB" "177MTCB" "474WTCB" "171HetCX" "174MTCX" "175HetCX" "177MTCX" "474WTCX")
  

for x in "${x[@]}"; do
    Rscript scripts/RPKM_TPM_featurecounts.R output/featurecounts/${x}.txt output/tpm/${x}
done
