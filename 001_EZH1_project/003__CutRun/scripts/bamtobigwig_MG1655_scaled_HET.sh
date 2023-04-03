#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.07929383602633

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.92654099341712

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.18237582286056

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.090065828845

bamCoverage --bam output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.0094484954298

bamCoverage --bam output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.39088014788949


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.15040566909726

bamCoverage --bam output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_HET_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

