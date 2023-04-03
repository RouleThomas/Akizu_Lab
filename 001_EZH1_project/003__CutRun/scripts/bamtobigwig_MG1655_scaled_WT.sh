#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.99640933572711

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 5.95347097546379

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.02154398563734

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 6.26256732495512


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.26091198521105

bamCoverage --bam output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 7.21351545650611


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.8490294751977


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_WT_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 9.98243812262504


