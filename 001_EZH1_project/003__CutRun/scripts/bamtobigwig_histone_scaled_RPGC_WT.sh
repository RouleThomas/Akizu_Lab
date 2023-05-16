#!/bin/bash
#SBATCH --mem=250G




bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.69085051546392 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.92551546391753 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.05154639175258 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398



bamCoverage --bam output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.09368556701031 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.1870882740448 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398



bamCoverage --bam output/bowtie2/8wN_WT_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.80237154150198 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398

bamCoverage --bam output/bowtie2/8wN_WT_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.20026350461133 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398


bamCoverage --bam output/bowtie2/8wN_WT_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_RPGC/8wN_WT_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 5.21607378129117 \
    --normalizeUsing RPGC \
    --effectiveGenomeSize 2913022398

