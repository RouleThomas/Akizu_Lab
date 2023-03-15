#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=72:00:00

module load sam-bcf-tools



x=("2dN_WT_R1" "2dN_WT_R2" "2dN_WT_R3"
   "2dN_KO_R1" "2dN_KO_R2" "2dN_KO_R3"
   "2dN_HET_R1" "2dN_HET_R2" "2dN_HET_R3")


for x in "${x[@]}"; do
    samtools index output/STAR/raw/${x}_Aligned.sortedByCoord.out.bam
    samtools index output/STAR/fastp/${x}_Aligned.sortedByCoord.out.bam
done
