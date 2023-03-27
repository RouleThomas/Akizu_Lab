#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL

x=("8wN_iPSCpatient_IGG_R1" "8wN_iPSCpatient_H3K27me3_R1"
   "8wN_iPSCpatient_IGG_R2" "8wN_iPSCpatient_H3K27me3_R2")


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/raw input/${x}_1.fq.gz
    fastqc -o output/fastqc/raw input/${x}_2.fq.gz
done