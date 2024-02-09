#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
"NPC_KOEF1aEZH1_IGG"
"NPC_KOEF1aEZH1_SUZ12"
"NPC_WT_EZH1cs"
"NPC_WT_EZH2"
"NPC_WT_H3K27ac"
"NPC_WT_H3K4me3"
"NPC_WT_SUZ12"
)


   

for x in "${x[@]}"; do
    fastqc -o output/fastqc/fastp output/fastp/${x}_1.fq.gz
    fastqc -o output/fastqc/fastp output/fastp/${x}_2.fq.gz
done