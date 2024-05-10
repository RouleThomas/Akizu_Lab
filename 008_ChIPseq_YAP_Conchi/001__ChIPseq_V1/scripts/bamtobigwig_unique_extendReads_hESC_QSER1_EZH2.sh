#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00


bamCoverage --bam output/bowtie2/hESC_WT_QSER1_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_QSER1_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 134 \
    --scaleFactor 1.410107084


bamCoverage --bam output/bowtie2/hESC_WT_QSER1_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_QSER1_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 132 \
    --scaleFactor 1.017742716


bamCoverage --bam output/bowtie2/hESC_YAPKO_QSER1_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_YAPKO_QSER1_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 132 \
    --scaleFactor 1.166519627



bamCoverage --bam output/bowtie2/hESC_YAPKO_QSER1_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_YAPKO_QSER1_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 132 \
    --scaleFactor 0.698720712





bamCoverage --bam output/bowtie2/hESC_WT_EZH2_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_EZH2_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 137 \
    --scaleFactor 1.354656714




bamCoverage --bam output/bowtie2/hESC_WT_EZH2_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_EZH2_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 152 \
    --scaleFactor 1.065553036


bamCoverage --bam output/bowtie2/hESC_YAPKO_EZH2_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_YAPKO_EZH2_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 129 \
    --scaleFactor 0.79313903




bamCoverage --bam output/bowtie2/hESC_YAPKO_EZH2_R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_YAPKO_EZH2_R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 132 \
    --scaleFactor 0.934339043





bamCoverage --bam output/bowtie2/hESC_WT_input_R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DiffBindTMM/hESC_WT_input_R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads 128 \
    --scaleFactor 1






