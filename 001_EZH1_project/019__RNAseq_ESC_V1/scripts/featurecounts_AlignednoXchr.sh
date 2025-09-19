#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


x=(
    "ESC_WT_R1"
    "ESC_KO_R1"
    "ESC_OEKO_R1"
    "ESC_WT_R2"
    "ESC_KO_R2"
    "ESC_OEKO_R2"
    "ESC_WT_R3"
    "ESC_KO_R3"
    "ESC_OEKO_R3"
    )
        
for x in "${x[@]}"; do
featureCounts -p -C -O -M --fraction -s 2 \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/${x}.txt output/STAR/fastp/${x}_AlignednoXchr.sortedByCoord.out.bam
done

