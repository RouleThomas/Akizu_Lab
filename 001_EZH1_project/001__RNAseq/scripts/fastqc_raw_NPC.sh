#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL



x=("NPC_WT_R1" "NPC_WT_R2" "NPC_WT_R3"
   "NPC_KO_R1" "NPC_KO_R2" "NPC_KO_R3"
   "NPC_HET_R1" "NPC_HET_R2" "NPC_HET_R3")
        
for x in "${x[@]}"; do
    fastqc -o output/fastqc input/${x}_1.fq.gz
    fastqc -o output/fastqc input/${x}_2.fq.gz
done