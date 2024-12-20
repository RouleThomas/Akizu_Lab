#!/bin/bash
#SBATCH --mem=500G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_WT_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.3135687

bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_WT_H3K27me3_010R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.1511360

bamCoverage --bam output/bowtie2/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_WT_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8573691


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_KO_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.4493761


bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_KO_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7184051

bamCoverage --bam output/bowtie2/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_KO_H3K27me3_014R2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.3219425


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_KOEF1aEZH1_H3K27me3_005R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.9534481


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_KOEF1aEZH1_H3K27me3_006R.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.2144370

bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_DBA_NORM_TMM/PSC_KOEF1aEZH1_H3K27me3_013R1.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8873456


