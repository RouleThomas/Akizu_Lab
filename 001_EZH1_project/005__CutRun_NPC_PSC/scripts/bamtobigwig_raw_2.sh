#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


input_list=("NPC_WT_EZH2"
"NPC_WT_H3K27me1"
"NPC_WT_H3K27me3"
"NPC_WT_H3K4me3"
"NPC_WT_IGG"
"NPC_WT_SUZ12"
"PSC_KOEF1aEZH1_EZH1cs"
"PSC_KOEF1aEZH1_EZH1pt"
"PSC_KOEF1aEZH1_H3K27me3"
"PSC_KOEF1aEZH1_HA")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

