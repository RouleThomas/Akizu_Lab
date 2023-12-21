#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



bamCoverage --bam output/bowtie2/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.444312079


bamCoverage --bam output/bowtie2/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.370942789


bamCoverage --bam output/bowtie2/50dN_WTQ731E_H3K27me3_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/50dN_WTQ731E_H3K27me3_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.386006703



