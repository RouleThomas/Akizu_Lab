#!/bin/bash
#SBATCH --mem=20G
#SBATCH --time=72:00:00


x=("HP14_Het" "HP42_Het" "HP43_Het" "HP20_KO" "HP38_KO" "HP41_KO"
   "CT14_Het" "CT42_Het" "CT43_Het" "CT20_KO" "CT38_KO" "CT41_KO"
   "CB14_Het" "CB42_Het" "CB43_Het" "CB20_KO" "CB38_KO" "CB41_KO")
  

for x in "${x[@]}"; do
    Rscript scripts/RPKM_TPM_featurecounts.R output/featurecounts/${x}.txt output/tpm/${x}
done
