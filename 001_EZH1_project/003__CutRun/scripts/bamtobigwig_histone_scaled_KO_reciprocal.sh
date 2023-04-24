#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G




bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7653614754906817


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.2624991543197346


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8857436365711714



bamCoverage --bam output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7263157894736834


bamCoverage --bam output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7014787430683908

bamCoverage --bam output/bowtie2/8wN_KO_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_KO_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.9730769230769262

