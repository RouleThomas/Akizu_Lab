#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00



bamCoverage --bam output/bowtie2/50dN_KOEF1aEZH1_H3K27me3_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.093518738

bamCoverage --bam output/bowtie2/50dN_KOEF1aEZH1_H3K27me3_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/50dN_KOEF1aEZH1_H3K27me3_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.636181604





