#!/bin/bash
#SBATCH --mail-user=roule@upenn.edu
#SBATCH --mail-type=ALL
#SBATCH --mem=50G
#SBATCH --time=100:00:00

module load R/4.2.2

bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R1.dupmark.sorted.fragments.bedgraph output/bowtie2/8wN_HET_IGG_R1.dupmark.sorted.fragments.bedgraph norm stringent output/SEACR/8wN_HET_H3K27me3_R1
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R2.dupmark.sorted.fragments.bedgraph output/bowtie2/8wN_HET_IGG_R2.dupmark.sorted.fragments.bedgraph norm stringent output/SEACR/8wN_HET_H3K27me3_R2
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R3.dupmark.sorted.fragments.bedgraph output/bowtie2/8wN_HET_IGG_R3.dupmark.sorted.fragments.bedgraph norm stringent output/SEACR/8wN_HET_H3K27me3_R3
bash scripts/SEACR_1.3.sh output/bowtie2/8wN_HET_H3K27me3_R4.dupmark.sorted.fragments.bedgraph output/bowtie2/8wN_HET_IGG_R4.dupmark.sorted.fragments.bedgraph norm stringent output/SEACR/8wN_HET_H3K27me3_R4
   
  

