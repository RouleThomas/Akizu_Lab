#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=72:00:00


x=("CB14_Het" "CB42_Het" "CB43_Het" "CB20_KO" "CB38_KO" "CB41_KO")
        
for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta_mice/ENCFF871VGR.gtf \
	-o output/featurecounts/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done