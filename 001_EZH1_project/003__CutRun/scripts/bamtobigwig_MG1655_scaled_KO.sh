#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.22396768402154

bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.44389587073609


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 7.20436864153202

bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.77749820273185

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.99388928828181

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.34148094895758


bamCoverage --bam output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_KO_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.22712334394577

