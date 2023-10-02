#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



x=(
"NPC_WT_EZH2"
"NPC_WT_H3K27me1"
"NPC_WT_H3K27me3"
"NPC_WT_H3K4me3"
"NPC_WT_IGG"
"NPC_WT_SUZ12"
"PSC_KOEF1aEZH1_EZH1cs"
"PSC_KOEF1aEZH1_EZH1pt"
"PSC_KOEF1aEZH1_H3K27me3"
"PSC_KOEF1aEZH1_HA"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




