#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam WT_p35_CB_Rep1/outs/possorted_genome_bam.bam \
    --outFileName output/bigwig/WT_p35_CB_Rep1.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --scaleFactor 1


bamCoverage --bam WT_p35_CB_Rep2/outs/possorted_genome_bam.bam \
    --outFileName output/bigwig/WT_p35_CB_Rep2.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --scaleFactor 1


bamCoverage --bam WT_p35_CB_Rep3/outs/possorted_genome_bam.bam \
    --outFileName output/bigwig/WT_p35_CB_Rep3.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --scaleFactor 1


