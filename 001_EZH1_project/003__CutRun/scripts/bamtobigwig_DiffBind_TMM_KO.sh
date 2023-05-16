#!/bin/bash
#SBATCH --mem=350G




bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.502057669


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.487749425


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.449980991 


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.883122496 



bamCoverage --bam output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.502057669 


bamCoverage --bam output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.487749425

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.449980991 

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_KO_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.883122496 

