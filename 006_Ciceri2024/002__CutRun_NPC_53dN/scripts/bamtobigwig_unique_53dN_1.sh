#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
    "53dN_WT_H3K27me3_R1"
    "53dN_WT_H3K27me3_R2"
    "53dN_WT_H3K4me3_R1"
    "53dN_WT_H3K4me3_R2"
    "53dN_WT_H3K9me3_R1"
    "53dN_WT_H3K9me3_R2"
    "53dN_WT_IGG_R1"
    "53dN_WT_IGG_R2"
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

