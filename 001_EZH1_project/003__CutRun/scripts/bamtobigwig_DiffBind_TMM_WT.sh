#!/bin/bash
#SBATCH --mem=250G




bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.412253844 


bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.186123259

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.337020479 



bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.40842928 


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.412253844



bamCoverage --bam output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.186123259

bamCoverage --bam output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.337020479


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBind_TMM/8wN_WT_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.40842928

