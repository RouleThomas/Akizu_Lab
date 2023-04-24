#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G




bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.4408840406795075

bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.6776700724827513


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.6226415094339607


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8470982142857131