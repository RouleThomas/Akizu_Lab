#!/bin/bash
#SBATCH --mem=500G




bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398



