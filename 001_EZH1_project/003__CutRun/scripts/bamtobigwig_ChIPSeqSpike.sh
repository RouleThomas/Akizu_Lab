#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=350G
#SBATCH --time=100:00:00


input_list=("8wN_HET_IGG_R1" "8wN_HET_H3K27me3_R1"
   "8wN_HET_IGG_R2" "8wN_HET_H3K27me3_R2"
   "8wN_HET_IGG_R3" "8wN_HET_H3K27me3_R3"
   "8wN_HET_IGG_R4" "8wN_HET_H3K27me3_R4"
   "8wN_KO_IGG_R1" "8wN_KO_H3K27me3_R1"
   "8wN_KO_IGG_R2" "8wN_KO_H3K27me3_R2"
   "8wN_KO_IGG_R3" "8wN_KO_H3K27me3_R3"
   "8wN_KO_IGG_R4" "8wN_KO_H3K27me3_R4"
   "8wN_WT_IGG_R1" "8wN_WT_H3K27me3_R1"
   "8wN_WT_IGG_R2" "8wN_WT_H3K27me3_R2" 
   "8wN_WT_IGG_R3" "8wN_WT_H3K27me3_R3"
   "8wN_WT_IGG_R4" "8wN_WT_H3K27me3_R4"
   "8wN_iPSCpatient_IGG_R1" "8wN_iPSCpatient_H3K27me3_R1"
   "8wN_iPSCpatient_IGG_R2" "8wN_iPSCpatient_H3K27me3_R2")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.dupmark.sorted.bam \
        --outFileName output/ChIPSeqSpike/${x}.bw \
        --outFileFormat bigwig \
        --binSize 50 \
        --numberOfProcessors 7 \
        --extendReads \
        --normalizeUsing RPGC
         --effectiveGenomeSize 2913022398
done




