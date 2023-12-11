#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"PSC_WT_EZH1cs_01FA"
"PSC_WT_EZH1cs_1FA"
"PSC_WT_H3K27me1_01FA"
"PSC_WT_H3K27me1_1FA"
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

