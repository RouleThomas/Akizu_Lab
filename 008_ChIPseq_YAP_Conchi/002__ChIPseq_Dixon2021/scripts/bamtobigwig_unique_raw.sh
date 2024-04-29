#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6   # Number of CPU cores per task



input_list=(
    "hESC_WT_DNMT3A_R1"
    "hESC_WT_DNMT3A_R2"
    "hESC_WT_DNMT3B_R1"
    "hESC_WT_DNMT3B_R2"
    "hESC_WT_H3K27me3_R1"
    "hESC_WT_H3K27me3_R2"
    "hESC_WT_H3K4me3_R1"
    "hESC_WT_H3K4me3_R2"
    "hESC_WT_input_R1"
    "hESC_WT_input_R2"
)

for x in "${input_list[@]}"; do
    bamCoverage --bam output/bowtie2/${x}.unique.dupmark.sorted.bam \
        --outFileName output/bigwig/${x}.unique.dupmark.sorted.bw \
        --outFileFormat bigwig \
        --binSize 1 \
        --numberOfProcessors 5 \
        --extendReads \
        --scaleFactor 0.5
done

