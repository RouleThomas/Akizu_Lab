#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=("8wN_KO_IGG_R1" "8wN_KO_H3K27me3_R1"
   "8wN_KO_IGG_R2" "8wN_KO_H3K27me3_R2"
   "8wN_KO_IGG_R3" "8wN_KO_H3K27me3_R3"
   "8wN_KO_IGG_R4" "8wN_KO_H3K27me3_R4")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.pgbam \
        --outFileName output/bigwig/${x}.pgbw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

