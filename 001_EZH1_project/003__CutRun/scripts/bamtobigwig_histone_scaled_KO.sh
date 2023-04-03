#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G




bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.30657216494845


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.80953608247423


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.12899484536082



bamCoverage --bam output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.3768115942029


bamCoverage --bam output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.42555994729908

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.02766798418972

