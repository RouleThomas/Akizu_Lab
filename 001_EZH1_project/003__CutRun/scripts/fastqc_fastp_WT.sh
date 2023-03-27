#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

x=("8wN_WT_IGG_R1" "8wN_WT_H3K27me3_R1"
   "8wN_WT_IGG_R2" "8wN_WT_H3K27me3_R2" 
   "8wN_WT_IGG_R3" "8wN_WT_H3K27me3_R3"
   "8wN_WT_IGG_R4" "8wN_WT_H3K27me3_R4")


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done