#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=72:00:00


x=("CT14_Het" "CT42_Het" "CT43_Het" "CT20_KO" "CT38_KO" "CT41_KO")
        
for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta_mice/ENCFF871VGR.gtf \
	-o output/featurecounts/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done