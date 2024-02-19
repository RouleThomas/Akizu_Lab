#!/bin/bash
#SBATCH --mem=150G
#SBATCH --time=200:00:00


input_list=("PSC_WT_EZH2"
"PSC_WT_IGG"
"PSC_WT_H3K27me3"
"PSC_WT_EZH1cs"
"PSC_WT_HA"
"PSC_WT_SUZ12")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.bam \
        --outFileName output/bigwig/${x}.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

