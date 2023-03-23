#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL



x=("NPC_HET_H3K27me3_R1"
"NPC_HET_H3K27me3_R2"
"NPC_HET_input_R1"
"NPC_HET_input_R2"
"NPC_KO_H3K27me3_R1"
"NPC_KO_H3K27me3_R2"
"NPC_KO_input_R1"
"NPC_KO_input_R2"
"NPC_WT_H3K27me3_R1"
"NPC_WT_H3K27me3_R2"
"NPC_WT_input_R1"
"NPC_WT_input_R2")


for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done