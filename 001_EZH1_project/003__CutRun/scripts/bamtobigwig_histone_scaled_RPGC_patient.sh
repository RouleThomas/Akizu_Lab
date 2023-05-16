#!/bin/bash
#SBATCH --mem=250G




bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.26817010309278 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.47564432989691 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.60606060606061 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.18050065876153 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398