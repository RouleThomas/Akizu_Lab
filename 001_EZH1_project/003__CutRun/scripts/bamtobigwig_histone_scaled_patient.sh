#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G




bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.26817010309278


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.47564432989691


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.60606060606061


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.18050065876153