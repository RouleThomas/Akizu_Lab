#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7




input_list=(
    "NPC_WT_H3K4me3_008"
    "NPC_WT_IGG_008"
    "NPC_WT_H3K4me3_005"
    "NPC_WT_IGG_005"
    "NPC_KO_H3K4me3_008"
    "NPC_KO_IGG_008"
    "NPC_KO_H3K4me3_005"
    "NPC_KO_IGG_005"
    "NPC_WT_H3K27me3_008"
    "NPC_WT_IGG_008"
    "NPC_WT_H3K27me3_005"
    "NPC_WT_IGG_005"
    "NPC_KO_H3K27me3_008"
    "NPC_KO_IGG_008"
    "NPC_KO_H3K27me3_005"
    "NPC_KO_IGG_005"
)





for x in "${input_list[@]}"; do
    bigWigToBedGraph output/bigwig/${x}.unique.dupmark.sorted.bw output/bigwig/${x}.unique.dupmark.sorted.bedGraph
done







