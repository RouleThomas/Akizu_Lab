#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/Ser5P_RNAPII.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/Ser5P_RNAPII.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 157 \
    --scaleFactor 0.5


