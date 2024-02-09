#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"50dNFA_KOEF1aEZH1_EZH1cs"
"50dNnative_KOEF1aEZH1_EZH1cs"
"50dNFA_KOEF1aEZH1_EZH2"
"50dNnative_KOEF1aEZH1_EZH2"
"50dNFA_KOEF1aEZH1_H3K27me3"
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

