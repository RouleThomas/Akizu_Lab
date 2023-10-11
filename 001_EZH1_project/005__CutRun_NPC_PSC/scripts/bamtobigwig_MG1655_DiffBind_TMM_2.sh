#!/bin/bash
#SBATCH --mem=250G




bamCoverage --bam output/bowtie2/NPC_WT_EZH2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_EZH2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.266595568




bamCoverage --bam output/bowtie2/NPC_WT_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.170001033





bamCoverage --bam output/bowtie2/NPC_WT_H3K4me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K4me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.21270728





bamCoverage --bam output/bowtie2/NPC_WT_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1



bamCoverage --bam output/bowtie2/NPC_WT_SUZ12.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_SUZ12.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.514593869



bamCoverage --bam output/bowtie2/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/PSC_KOEF1aEZH1_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.445496541

