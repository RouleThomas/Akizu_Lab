#!/bin/bash
#SBATCH --mem=150G



input_list=("PSC_KOEF1aEZH1_IGG"
"PSC_KOEF1aEZH1_SUZ12"
"PSC_KOsynEZH1_EZH1cs"
"PSC_KOsynEZH1_EZH1pt"
"PSC_KOsynEZH1_H3K27me3"
"PSC_KOsynEZH1_HA"
"PSC_KOsynEZH1_IGG"
"PSC_KOsynEZH1_SUZ12")

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads \
        --scaleFactor 0.5
done

