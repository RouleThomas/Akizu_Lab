#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



bamCoverage --bam output/bowtie2/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.052366383

bamCoverage --bam output/bowtie2/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.04068647


bamCoverage --bam output/bowtie2/50dN_WTQ731E_H3K27me3_R3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/50dN_WTQ731E_H3K27me3_R3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.995054876



