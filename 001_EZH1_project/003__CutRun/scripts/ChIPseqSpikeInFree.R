#! /usr/bin/env Rscript


# Load packages
library("Rsamtools")
library("GenomicAlignments")
library("ChIPseqSpikeInFree")

# Create sample_meta.txt; tab delimited format `output/ChIPseqSpikeInFree/sample_meta.txt
metaFile <- "output/ChIPseqSpikeInFree/sample_meta.txt"

bams <- c("output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.bam","output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.bam","output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.bam","output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.bam","output/bowtie2/8wN_iPSCpatient_H3K27me3_R1.dupmark.sorted.bam","output/bowtie2/8wN_KO_H3K27me3_R1.dupmark.sorted.bam","output/bowtie2/8wN_KO_H3K27me3_R2.dupmark.sorted.bam","output/bowtie2/8wN_KO_H3K27me3_R3.dupmark.sorted.bam","output/bowtie2/8wN_KO_H3K27me3_R4.dupmark.sorted.bam","output/bowtie2/8wN_WT_H3K27me3_R1.dupmark.sorted.bam","output/bowtie2/8wN_WT_H3K27me3_R2.dupmark.sorted.bam","output/bowtie2/8wN_WT_H3K27me3_R3.dupmark.sorted.bam","output/bowtie2/8wN_WT_H3K27me3_R4.dupmark.sorted.bam")

# Run ChIPSpikeInFree
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "hg38", metaFile = metaFile, prefix = "H3K27me3")