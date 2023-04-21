#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=75G




bamCoverage --bam output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.760228354


bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.663115954


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.926914153


bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.557530594


bamCoverage --bam output/bowtie2/8wN_KO_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.726315789


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.765361475



bamCoverage --bam output/bowtie2/8wN_KO_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_KO_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.735028712




bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.650588035


bamCoverage --bam output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_H3K27me3_R3.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1



bamCoverage --bam output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.bam \
    --outFileName output/bigwig_histone/8wN_HET_IGG_R4.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.411219763


