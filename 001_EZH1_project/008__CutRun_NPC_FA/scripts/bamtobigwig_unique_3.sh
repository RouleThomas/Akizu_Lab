#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"NPC_KOEF1aEZH1_IGG"
"NPC_KOEF1aEZH1_SUZ12"
"NPC_WT_EZH1cs"
"NPC_WT_EZH2"
"NPC_WT_H3K27ac"
"NPC_WT_H3K4me3"
"NPC_WT_SUZ12"
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

