#!/bin/bash
#SBATCH --mem=250G




bamCoverage --bam output/bowtie2/NPC_KO_EZH2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_EZH2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.792209099



bamCoverage --bam output/bowtie2/NPC_KO_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.870425955


bamCoverage --bam output/bowtie2/NPC_KO_H3K4me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.830386499


bamCoverage --bam output/bowtie2/NPC_KO_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

bamCoverage --bam output/bowtie2/NPC_KO_SUZ12.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_SUZ12.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.713876255






























