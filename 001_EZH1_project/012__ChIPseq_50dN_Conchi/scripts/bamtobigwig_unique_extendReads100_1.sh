#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
"50dN_KO_input_R1"
"50dN_KOEF1aEZH1_input_R1"
"50dN_WT_input_R1"
"50dN_WT_EZH1_R1"
"50dN_WT_EZH1_R2"
"50dN_KO_EZH2_R1"
"50dN_KO_EZH2_R2"
"50dN_KOEF1aEZH1_EZH2_R1"
"50dN_KOEF1aEZH1_EZH2_R2"
"50dN_WT_EZH2_R1"
"50dN_WT_EZH2_R2"
"50dN_WT_H3K27me1_R1"
"50dN_WT_H3K27me1_R2"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig_extendReads100/${x}.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 7 \
        --extendReads 100 \
        --scaleFactor 0.5
done

