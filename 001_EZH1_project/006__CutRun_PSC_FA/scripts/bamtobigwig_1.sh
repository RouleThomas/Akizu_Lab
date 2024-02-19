#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



input_list=("PSC_KOEF1aEZH1_H3K27me3"
"PSC_KOEF1aEZH1_HA"
"PSC_KOEF1aEZH1_EZH1cs"
"PSC_KOEF1aEZH1_EZH2"
"PSC_KOEF1aEZH1_SUZ12"
"PSC_KOEF1aEZH1_IGG")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

