#!/bin/bash
#SBATCH --mem=250G




bamCoverage --bam output/bowtie2/NPC_WT_EZH2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_WT_EZH2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.30168537



bamCoverage --bam output/bowtie2/NPC_WT_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_WT_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.167456178




bamCoverage --bam output/bowtie2/NPC_WT_H3K4me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_WT_H3K4me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.2446728




bamCoverage --bam output/bowtie2/NPC_WT_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_WT_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1



bamCoverage --bam output/bowtie2/NPC_WT_SUZ12.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_WT_SUZ12.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.574337542


bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.315573789
