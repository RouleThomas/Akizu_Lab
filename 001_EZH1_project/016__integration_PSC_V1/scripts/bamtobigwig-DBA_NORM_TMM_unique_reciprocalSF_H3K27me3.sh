#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.761284887


bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.868707086

bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.166358806


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 2.225307487


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.391972301

bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.756462554


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.048824787


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.823426822

bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM_reciprocalSF/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.126956622


