#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"PSC_KOEF1aEZH1_IGG"
"PSC_KOEF1aEZH1_SUZ12"
"PSC_KOsynEZH1_EZH1cs"
"PSC_KOsynEZH1_EZH1pt"
"PSC_KOsynEZH1_H3K27me3"
"PSC_KOsynEZH1_HA"
"PSC_KOsynEZH1_IGG"
"PSC_KOsynEZH1_SUZ12"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done