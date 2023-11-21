#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=72:00:00


x=("S_CB_KO1" "S_CB_KO2" "S_CB_KO3" "S_CB_WT1" "S_CB_WT2" "S_CB_WT3")
        
for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta_mice/ENCFF871VGR.gtf \
	-o output/featurecounts/${x}.txt output/bam/1month_rnaseqyz072420/${x}.bam
done