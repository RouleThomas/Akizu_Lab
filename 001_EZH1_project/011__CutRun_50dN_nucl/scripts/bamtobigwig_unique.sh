#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"50dN_WT_EZH1"
"50dN_WT_EZH2"
"50dN_WT_H3K27ac"
"50dN_WT_H3K27me1AM"
"50dN_WT_H3K27me1OR"
"50dN_WT_H3K27me3"
"50dN_WT_IGG"
"50dN_WT_SUZ12"
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

