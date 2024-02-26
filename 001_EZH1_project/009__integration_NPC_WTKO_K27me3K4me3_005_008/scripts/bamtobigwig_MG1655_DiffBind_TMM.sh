#!/bin/bash
#SBATCH --mem=250G
#SBATCH --time=200:00:00
#SBATCH --cpus-per-task=7



bamCoverage --bam output/bowtie2/NPC_WT_IGG_008.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_IGG_008.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/NPC_WT_IGG_005.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_IGG_005.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1

bamCoverage --bam output/bowtie2/NPC_WT_H3K4me3_008.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K4me3_008.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7577088

bamCoverage --bam output/bowtie2/NPC_WT_H3K4me3_005.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K4me3_005.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.0857384

bamCoverage --bam output/bowtie2/NPC_WT_H3K27me3_008.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K27me3_008.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.1712504



bamCoverage --bam output/bowtie2/NPC_WT_H3K27me3_005.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_WT_H3K27me3_005.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.1253542

bamCoverage --bam output/bowtie2/NPC_KO_IGG_008.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_IGG_008.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/NPC_KO_IGG_005.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_IGG_005.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1


bamCoverage --bam output/bowtie2/NPC_KO_H3K4me3_008.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3_008.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.7061478


bamCoverage --bam output/bowtie2/NPC_KO_H3K4me3_005.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K4me3_005.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.7403885


bamCoverage --bam output/bowtie2/NPC_KO_H3K27me3_008.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27me3_008.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 1.0903367

bamCoverage --bam output/bowtie2/NPC_KO_H3K27me3_005.unique.dupmark.sorted.bam \
    --outFileName output/bigwig_MG1655_DiffBind_TMM/NPC_KO_H3K27me3_005.unique.dupmark.sorted.bw \
    --outFileFormat bigwig \
    --binSize 1 \
    --numberOfProcessors 7 \
    --extendReads \
    --scaleFactor 0.8471138




