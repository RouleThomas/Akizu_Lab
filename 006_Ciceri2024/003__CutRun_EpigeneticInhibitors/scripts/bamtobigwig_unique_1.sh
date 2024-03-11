#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00



input_list=(
    "EZH2inh_H3K27me3_R1"
    "EHMTinh_noAB_R2"
    "EHMTinh_noAB_R1"
    "EHMTinh_H3K4me3_R2"
    "EHMTinh_H3K4me3_R1"
    "EHMTinh_H3K27me3_R2"
    "EHMTinh_H3K27me3_R1"
    "DOT1Linh_noAB_R2"
    "DOT1Linh_noAB_R1"
    "DOT1Linh_H3K4me3_R2"
    "DOT1Linh_H3K4me3_R1"
    "DOT1Linh_H3K27me3_R2"
    "DOT1Linh_H3K27me3_R1"
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

