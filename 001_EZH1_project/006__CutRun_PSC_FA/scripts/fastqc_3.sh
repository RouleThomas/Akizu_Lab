#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"PSC_WT_EZH2"
"PSC_WT_IGG"
"PSC_WT_H3K27me3"
"PSC_WT_EZH1cs"
"PSC_WT_HA"
"PSC_WT_SUZ12"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/raw input/${x}_1.fq.gz
    fastqc -o output/fastqc/raw input/${x}_2.fq.gz
done