#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G



bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_iPSCpatient_H3K27me3_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.0050867743866

bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_iPSCpatient_IGG_R2.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.31067063777344

bamCoverage --bam output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 4.18940754039497


bamCoverage --bam output/bowtie2/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655/8wN_iPSCpatient_IGG_R1.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 10 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.18445106295574