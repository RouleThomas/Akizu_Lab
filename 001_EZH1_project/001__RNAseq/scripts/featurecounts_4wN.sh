#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=20G
#SBATCH --time=72:00:00


x=("4wN_WT_R1" "4wN_WT_R2" "4wN_KO_R1"
   "4wN_KO_R2" "4wN_HET_R1" "4wN_HET_R2"
   "4wN_HET_R3" "4wN_HET_R4" "4wN_iPSCWT_R1"
   "4wN_iPSCWT_R2" "4wN_iPSCpatient_R1" "4wN_iPSCpatient_R2")
        
for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/gencode.v19.annotation.gtf \
	-o output/featurecounts/${x}.txt output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done