#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"NPC_WT_IGG"
"NPC_WT_H3K27me3"
"NPC_KO_H3K4me3"
"NPC_KO_H3K27me3"
"NPC_KO_SUZ12"
"NPC_KOEF1aEZH1_H3K27ac"
"NPC_KOEF1aEZH1_EZH2"
"50dNFA_KOEF1aEZH1_IGG"
"50dNnative_KOEF1aEZH1_IGG"
"50dNnative_KOEF1aEZH1_H3K27me3"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done




