#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"50dN_KO_EZH1cs_R1"
"50dN_KO_EZH1cs_R2"
"50dN_KO_EZH2_R1"
"50dN_KO_EZH2_R2"
"50dN_KO_H3K27me3_R1"
"50dN_KO_H3K27me3_R2"
"50dN_KO_IGG_R1"
"50dN_KO_IGG_R2"
"50dN_KO_SUZ12_R1"
"50dN_KO_SUZ12_R2"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done