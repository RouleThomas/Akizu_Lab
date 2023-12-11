#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"50dN_KOEF1aEZH1_EZH1cs_R1"
"50dN_KOEF1aEZH1_EZH1cs_R2"
"50dN_KOEF1aEZH1_EZH2_R1"
"50dN_KOEF1aEZH1_EZH2_R2"
"50dN_KOEF1aEZH1_H3K27me3_R1"
"50dN_KOEF1aEZH1_H3K27me3_R2"
"50dN_KOEF1aEZH1_IGG_R1"
"50dN_KOEF1aEZH1_SUZ12_R1"
"50dN_KOEF1aEZH1_SUZ12_R2"
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

