#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G




bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.5914183370169948


bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.3418201039555981

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.4874371859296476



bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.3232390552755446


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.4572289156626497



bamCoverage --bam output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.3568406205923832

bamCoverage --bam output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8331503841931948


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_WT_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.191715079565547

