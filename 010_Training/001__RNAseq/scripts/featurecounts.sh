#!/bin/bash
#SBATCH --mem=25G
#SBATCH --time=72:00:00


x=("WT_Rep1" "WT_Rep2" "WT_Rep3" "KO_Rep1" "KO_Rep2" "KO_Rep3")

        
for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /Master/meta/gencode.v47.annotation.gtf \
	-o output/featurecounts/${x}.txt output/STAR/${x}_Aligned.sortedByCoord.out.bam
done



