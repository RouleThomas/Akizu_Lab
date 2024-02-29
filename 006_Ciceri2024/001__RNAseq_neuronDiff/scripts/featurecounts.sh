#!/bin/bash
#SBATCH --mem=50G
#SBATCH --time=100:00:00


x=(
    "100dN_WT_R1"
    "75dN_WT_R1"
    "50dN_WT_R1"
    "25dN_WT_R1"
    "NPC_WT_R1"
    "ESC_WT_R1"
    "100dN_WT_R2"
    "75dN_WT_R2"
    "50dN_WT_R2"
    "25dN_WT_R2"
    "NPC_WT_R2"
    "ESC_WT_R2"
    "100dN_WT_R3"
    "75dN_WT_R3"
    "50dN_WT_R3"
    "25dN_WT_R3"
    "NPC_WT_R3"
    "ESC_WT_R3"
)
        
for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf \
	-o output/featurecounts/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done