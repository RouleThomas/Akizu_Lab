#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.58479381443299

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.88492268041237


bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.05090206185567

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.85734536082474


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.38471673254282


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.13570487483531


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.05270092226614


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.55994729907773


