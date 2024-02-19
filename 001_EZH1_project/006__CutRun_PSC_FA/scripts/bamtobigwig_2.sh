#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00



input_list=("PSC_KO_H3K27me3"
"PSC_KO_EZH2"
"PSC_KO_SUZ12"
"PSC_KO_IGG"
"PSC_KO_EZH1cs"
"PSC_KO_HA")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

