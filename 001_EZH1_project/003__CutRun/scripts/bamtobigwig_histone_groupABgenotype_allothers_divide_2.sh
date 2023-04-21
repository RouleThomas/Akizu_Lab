#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=75G




bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.367789654


bamCoverage --bam output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.701478743

bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.262499154


bamCoverage --bam output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.973076923


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.885743637

bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

