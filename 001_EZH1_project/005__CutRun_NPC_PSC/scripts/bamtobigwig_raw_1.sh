#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00


input_list=("NPC_KO_EZH1cs"
"NPC_KO_EZH1pt"
"NPC_KO_EZH2"
"NPC_KO_H3K27me1"
"NPC_KO_H3K27me3"
"NPC_KO_H3K4me3"
"NPC_KO_IGG"
"NPC_KO_SUZ12"
"NPC_WT_EZH1cs"
"NPC_WT_EZH1pt")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

