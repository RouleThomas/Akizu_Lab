#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam RNA_Bap1KO/outs/possorted_genome_bam.bam \
    --outFileName output/bigwig/RNA_Bap1KO.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --scaleFactor 1



