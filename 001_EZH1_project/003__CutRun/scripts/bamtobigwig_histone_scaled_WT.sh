#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G




bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.69085051546392


bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.92551546391753

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.05154639175258



bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.09368556701031


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.1870882740448



bamCoverage --bam output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.80237154150198

bamCoverage --bam output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.20026350461133


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_WT_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 5.21607378129117

