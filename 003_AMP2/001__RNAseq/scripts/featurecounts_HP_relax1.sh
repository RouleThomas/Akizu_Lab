#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=72:00:00


x=("HP14_Het" "HP42_Het" "HP43_Het" "HP20_KO" "HP38_KO" "HP41_KO")
        
for x in "${x[@]}"; do
featureCounts -p -C -O -M --fraction \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta_mice/ENCFF871VGR.gtf \
	-o output/featurecounts_relax1/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done
