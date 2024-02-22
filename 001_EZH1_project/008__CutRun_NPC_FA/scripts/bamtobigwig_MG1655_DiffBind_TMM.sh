#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam output/bowtie2/NPC_WT_H3K27ac.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K27ac.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.062330996


bamCoverage --bam output/bowtie2/NPC_WT_H3K4me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K4me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.695790606


bamCoverage --bam output/bowtie2/NPC_WT_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.099062335

bamCoverage --bam output/bowtie2/NPC_KO_H3K27ac.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27ac.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.139698939

bamCoverage --bam output/bowtie2/NPC_KO_H3K4me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.459619623



bamCoverage --bam output/bowtie2/NPC_KO_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.996524124




bamCoverage --bam output/bowtie2/NPC_WT_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/NPC_KO_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1




