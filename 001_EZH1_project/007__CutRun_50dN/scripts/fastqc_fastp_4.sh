#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"PSC_WT_EZH1cs_01FA"
"PSC_WT_EZH1cs_1FA"
"PSC_WT_H3K27me1_01FA"
"PSC_WT_H3K27me1_1FA"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done