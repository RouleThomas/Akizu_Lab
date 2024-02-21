#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"NPC_WT_IGG"
"NPC_WT_H3K27me3"
"NPC_KO_H3K4me3"
"NPC_KO_H3K27me3"
"NPC_KO_SUZ12"
"NPC_KOEF1aEZH1_H3K27ac"
"NPC_KOEF1aEZH1_EZH2"
"50dNFA_KOEF1aEZH1_IGG"
"50dNnative_KOEF1aEZH1_IGG"
"50dNnative_KOEF1aEZH1_H3K27me3"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

