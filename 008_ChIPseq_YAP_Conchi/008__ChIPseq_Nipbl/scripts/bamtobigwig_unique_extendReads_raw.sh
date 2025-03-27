#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/Nipbl.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/Nipbl.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 138 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/input.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/input.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --scaleFactor 0.5


    
