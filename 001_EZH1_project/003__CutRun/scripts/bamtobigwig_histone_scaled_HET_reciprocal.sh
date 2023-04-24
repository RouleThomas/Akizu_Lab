#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.6309969100666774

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.5305257400697344


bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.9515634580012263

bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_H3K27me3_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.3499751950570516


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7221693625118932


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8805104408352665


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.9499374217772212


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_NotGenotypeGroup/8wN_HET_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 50 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.3906330416881118


