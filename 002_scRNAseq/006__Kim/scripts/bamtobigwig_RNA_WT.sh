#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam RNA_WT/outs/possorted_genome_bam.bam \
    --outFileName output/bigwig/RNA_WT.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --scaleFactor 1



