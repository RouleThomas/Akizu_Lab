#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=6   # Number of CPU cores per task



input_list=(
    "hESC_WT_EZH2_R1"
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

