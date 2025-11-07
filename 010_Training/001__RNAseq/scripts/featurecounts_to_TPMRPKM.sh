#!/bin/bash
#SBATCH --mem=25G
#SBATCH --time=72:00:00

module load R/4.2.3



x=("WT_Rep1" "WT_Rep2" "WT_Rep3" "KO_Rep1" "KO_Rep2" "KO_Rep3")

  

for x in "${x[@]}"; do
    Rscript scripts/RPKM_TPM_featurecounts.R output/featurecounts/${x}.txt output/tpm/${x}
done
