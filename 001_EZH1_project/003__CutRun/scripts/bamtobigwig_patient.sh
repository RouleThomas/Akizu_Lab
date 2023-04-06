#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=75G



input_list=("8wN_iPSCpatient_IGG_R1" "8wN_iPSCpatient_H3K27me3_R1"
   "8wN_iPSCpatient_IGG_R2" "8wN_iPSCpatient_H3K27me3_R2")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done


