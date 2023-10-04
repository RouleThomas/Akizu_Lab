#!/bin/bash
#SBATCH --mem=250G



input_list=(
"NPC_KO_IGG"
"NPC_KO_SUZ12"
"NPC_WT_EZH1cs"
"NPC_WT_EZH1pt"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

