#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8337495

bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.6454013

bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.0568582


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.6979207


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.0157988

bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.0308526


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7405806


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8654724

bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_LIB/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.1133659


