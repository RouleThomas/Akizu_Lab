#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/hESC_WT_QSER1FLAG_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_QSER1FLAG_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.593124384


bamCoverage --bam output/bowtie2/hESC_WT_QSER1FLAG_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_QSER1FLAG_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 3.479459014



bamCoverage --bam output/bowtie2/hESC_WT_H3K27me3_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_H3K27me3_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.751266429



bamCoverage --bam output/bowtie2/hESC_WT_H3K27me3_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_H3K27me3_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.520177776




bamCoverage --bam output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_input_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1



bamCoverage --bam output/bowtie2/hESC_WT_input_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_input_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1









