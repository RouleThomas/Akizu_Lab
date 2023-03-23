#!/bin/bash
#SBATCH --mem=20G
#SBATCH --time=72:00:00

module load R/4.2.2



x=("2dN_WT_R1" "2dN_WT_R2" "2dN_WT_R3"
   "2dN_KO_R1" "2dN_KO_R2" "2dN_KO_R3"
   "2dN_HET_R1" "2dN_HET_R2" "2dN_HET_R3"
   "ESC_WT_R1" "ESC_WT_R2" "ESC_WT_R3"
   "ESC_KO_R1" "ESC_KO_R2" "ESC_KO_R3"
   "ESC_HET_R1" "ESC_HET_R2" "ESC_HET_R3"
   "NPC_WT_R1" "NPC_WT_R2" "NPC_WT_R3"
   "NPC_KO_R1" "NPC_KO_R2" "NPC_KO_R3"
   "NPC_HET_R1" "NPC_HET_R2" "NPC_HET_R3"
   "4wN_WT_R1" "4wN_WT_R2" "4wN_KO_R1"
   "4wN_KO_R2" "4wN_HET_R1" "4wN_HET_R2"
   "4wN_HET_R3" "4wN_HET_R4" "4wN_iPSCWT_R1"
   "4wN_iPSCWT_R2" "4wN_iPSCpatient_R1" "4wN_iPSCpatient_R2"
   "8wN_WT_R1" "8wN_WT_R2" "8wN_WT_R3" "8wN_WT_R4" "8wN_KO_R1"
   "8wN_KO_R2" "8wN_KO_R3" "8wN_KO_R4" "8wN_HET_R1" "8wN_HET_R2"
   "8wN_HET_R3" "8wN_HET_R4" "8wN_iPSCpatient_R1" "8wN_iPSCpatient_R2")
  

for x in "${x[@]}"; do
    Rscript scripts/RPKM_TPM_featurecounts.R output/featurecounts/${x}.txt output/tpm/${x}
done
