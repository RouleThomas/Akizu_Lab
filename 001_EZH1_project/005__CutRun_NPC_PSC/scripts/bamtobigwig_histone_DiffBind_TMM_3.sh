#!/bin/bash
#SBATCH --mem=250G







bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_HA.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOEF1aEZH1_HA.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1




bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOEF1aEZH1_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_SUZ12.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOEF1aEZH1_SUZ12.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.32043567


bamCoverage --bam output/bowtie2/PSC_KOsynEZH1_EZH1cs.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOsynEZH1_EZH1cs.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/PSC_KOsynEZH1_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOsynEZH1_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.764069188


bamCoverage --bam output/bowtie2/PSC_KOsynEZH1_HA.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOsynEZH1_HA.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/PSC_KOsynEZH1_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOsynEZH1_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1



bamCoverage --bam output/bowtie2/PSC_KOsynEZH1_SUZ12.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOsynEZH1_SUZ12.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.810146928






