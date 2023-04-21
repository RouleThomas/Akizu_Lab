#!/bin/bash
#SBATCH --mem=25G
#SBATCH --time=72:00:00


x=("8wN_WT_R1" "8wN_WT_R2" "8wN_WT_R3" "8wN_WT_R4" "8wN_KO_R1"
   "8wN_KO_R2" "8wN_KO_R3" "8wN_KO_R4" "8wN_HET_R1" "8wN_HET_R2"
   "8wN_HET_R3" "8wN_HET_R4" "8wN_iPSCpatient_R1" "8wN_iPSCpatient_R2"
   "2dN_WT_R1" "2dN_WT_R2" "2dN_WT_R3"
   "2dN_KO_R1" "2dN_KO_R2" "2dN_KO_R3"
   "2dN_HET_R1" "2dN_HET_R2" "2dN_HET_R3")
        
for x in "${x[@]}"; do
featureCounts -p -C -O \
	-a /scr1/users/roulet/Akizu_Lab/Master/meta/ENCFF159KBI.gtf \
	-o output/featurecounts_hg38/${x}.txt output/STAR_hg38/${x}_Aligned.sortedByCoord.out.bam
done