#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=72:00:00


x=("171HetCX" "174MTCX" "175HetCX" "177MTCX" "474WTCX")
        

for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta_mice/ENCFF871VGR.gtf \
	-o output/featurecounts/${x}.txt output/bam/1year_AL1804271_R2_new_analysis/${x}.bam
done