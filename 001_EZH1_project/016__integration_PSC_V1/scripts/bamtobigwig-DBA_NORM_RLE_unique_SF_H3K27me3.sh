#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.4792448

bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.2771297

bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.9315771


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.4861315


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7737804

bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.4222858


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.0176293


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.3667435

bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_RLE/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.9730354


