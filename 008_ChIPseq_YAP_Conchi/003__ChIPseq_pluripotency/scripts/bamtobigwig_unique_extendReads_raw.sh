#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/hESC_WT_input_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/hESC_WT_TEAD4_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/hESC_WT_TEAD4_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 132 \
    --scaleFactor 0.5


bamCoverage --bam output/bowtie2/hESC_WT_YAP1_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig/hESC_WT_YAP1_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 125 \
    --scaleFactor 0.5


    
