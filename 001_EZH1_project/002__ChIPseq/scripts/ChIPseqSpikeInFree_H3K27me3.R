#! /usr/bin/env Rscript


# Load packages
library("Rsamtools")
library("GenomicAlignments")
library("ChIPseqSpikeInFree")

# Create sample_meta.txt; tab delimited format `output/ChIPseqSpikeInFree/sample_meta.txt
metaFile <- "output/ChIPseqSpikeInFree_H3K27me3/sample_meta.txt"

bams <- c("output/bowtie2_endtoend/2dN_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/2dN_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_HET_H3K27me3_R2.dupmark.sorted.bam",  "output/bowtie2_endtoend/ESC_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/ESC_WT_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_HET_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_KO_H3K27me3_R2.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_H3K27me3_R1.dupmark.sorted.bam", "output/bowtie2_endtoend/NPC_WT_H3K27me3_R2.dupmark.sorted.bam")

# Run ChIPSpikeInFree
ChIPseqSpikeInFree(bamFiles = bams, chromFile = "hg38", metaFile = metaFile, prefix = "H3K27me3")