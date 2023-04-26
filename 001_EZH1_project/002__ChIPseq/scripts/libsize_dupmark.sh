#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=100:00:00



module load sam-bcf-tools/1.6


 

samtools flagstat output/bowtie2_endtoend/2dN_HET_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/2dN_HET_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/2dN_KO_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/2dN_KO_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/2dN_WT_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/2dN_WT_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/ESC_HET_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/ESC_HET_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/ESC_KO_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/ESC_KO_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/ESC_WT_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/ESC_WT_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/NPC_HET_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/NPC_HET_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/NPC_KO_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/NPC_KO_H3K27me3_R2.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/NPC_WT_H3K27me3_R1.dupmark.sorted.bam 
samtools flagstat output/bowtie2_endtoend/NPC_WT_H3K27me3_R2.dupmark.sorted.bam 


