#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

x=("ESC_HET_R1" "ESC_HET_R2" "ESC_HET_R3"
   "NPC_WT_R1" "NPC_WT_R2" "NPC_WT_R3"
   "NPC_KO_R1" "NPC_KO_R2" "NPC_KO_R3"
   "NPC_HET_R1" "NPC_HET_R2" "NPC_HET_R3")
        
for x in "${x[@]}"; do
    fastqc -o output/fastqc output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc output/fastp/${x}_2.fq.gz
done