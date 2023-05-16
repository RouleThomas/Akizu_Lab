#!/bin/bash
#SBATCH --mem=250G



bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.58479381443299 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398



bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.88492268041237 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.05090206185567 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.85734536082474 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.38471673254282 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.13570487483531 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.05270092226614 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_HET_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.55994729907773 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


