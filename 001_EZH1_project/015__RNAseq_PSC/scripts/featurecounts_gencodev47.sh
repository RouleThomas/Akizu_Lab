#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
    "PSC_WT_R1"
    "PSC_WT_R2"
    "PSC_WT_R3"
    "PSC_KO_R1"
    "PSC_KO_R2"
    "PSC_KO_R3"
    "PSC_KOEF1aEZH1_R1"
    "PSC_KOEF1aEZH1_R2"
    "PSC_KOEF1aEZH1_R3"
    )
        
for x in "${x[@]}"; do
featureCounts -p -C -O -M --fraction -s 2 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts_v47/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done

