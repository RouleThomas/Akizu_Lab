#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=100:00:00




x=(
    "53dN_WT_H3K27ac_R1"
    "53dN_WT_H3K27ac_R2"
    "53dN_WT_H3K27me3_R1"
    "53dN_WT_H3K27me3_R2"
    "53dN_WT_H3K4me3_R1"
    "53dN_WT_H3K4me3_R2"
    "53dN_WT_H3K9me3_R1"
    "53dN_WT_H3K9me3_R2"
    "53dN_WT_IGG_R1"
    "53dN_WT_IGG_R2"
    "NPC_WT_H3K27ac_R1"
    "NPC_WT_H3K27ac_R2"
    "NPC_WT_H3K27me3_R1"
    "NPC_WT_H3K27me3_R2"
    "NPC_WT_H3K4me3_R1"
    "NPC_WT_H3K4me3_R2"
    "NPC_WT_H3K9me3_R1"
    "NPC_WT_H3K9me3_R2"
    "NPC_WT_IGG_R1"
    "NPC_WT_IGG_R2"
)

for x in "${x[@]}"; do
    fastp -i input/${x}_1.fq.gz -I input/${x}_2.fq.gz \
    -o output/fastp/${x}_1.fq.gz -O output/fastp/${x}_2.fq.gz \
    -j output/fastp/${x} -h output/fastp/${x}
done