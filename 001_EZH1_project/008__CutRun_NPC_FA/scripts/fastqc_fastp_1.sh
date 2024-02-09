#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"50dNFA_KOEF1aEZH1_EZH1cs"
"50dNnative_KOEF1aEZH1_EZH1cs"
"50dNFA_KOEF1aEZH1_EZH2"
"50dNnative_KOEF1aEZH1_EZH2"
"50dNFA_KOEF1aEZH1_H3K27me3"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done