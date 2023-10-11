#!/bin/bash
#SBATCH --mem=250G



bamCoverage --bam output/bowtie2/NPC_KO_EZH1cs.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_KO_EZH1cs.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1



bamCoverage --bam output/bowtie2/NPC_KO_EZH2.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_KO_EZH2.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.814156687


bamCoverage --bam output/bowtie2/NPC_KO_H3K27me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_KO_H3K27me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.868532645


bamCoverage --bam output/bowtie2/NPC_KO_H3K4me3.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_KO_H3K4me3.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.852010831


bamCoverage --bam output/bowtie2/NPC_KO_IGG.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_KO_IGG.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

bamCoverage --bam output/bowtie2/NPC_KO_SUZ12.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_KO_SUZ12.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.74203529



bamCoverage --bam output/bowtie2/NPC_WT_EZH1cs.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_histone_DiffBind_TMM/NPC_WT_EZH1cs.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1






























